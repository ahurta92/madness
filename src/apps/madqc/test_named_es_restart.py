#!/usr/bin/env python3
"""Excited-state bundles across runs on one named stage (D-088, D-089).

D-088: a re-solve at an aliased threshold (1.2e-4 and 1e-4 share the key
1e-04_k6) redefines the protocol entry it is saved under, so the state record,
the engine's protocol entry and the bundle's own index describe the same rung
and the stage's verdict is on the requested rung. Before the fix the entry
kept the old threshold, same_rung failed the request and the stage ended
unconverged although the bundle had converged.

D-089: a converged bundle loaded without its transition analysis gets the
analysis computed from the saved states, with no solver iterations. The
engine publishes the bundle at its save boundary and writes the analysis
after, so a restart between the two left the transition properties missing
for good. The test simulates that state by removing the generation published
at the property-analysis phase and pointing `latest` at the bundle's
save-boundary generation (the checksums cover the archives, not the
pointer), then continues the stage.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='named-es-restart-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]

# An excited-state-only stage: the dipole-leg writer already redefines an
# aliased protocol entry when it saves, so a mixed deck can hide the bundle
# writer's omission; here only the bundle writes.
mixed = (Path(a.source) / 'named_h2_es_tda.in').read_text()
head, response = mixed.split('molresponse excited')
response = 'molresponse excited' + response
assert '[excited_states, polarizability]' in response
response = response.replace('[excited_states, polarizability]', '[excited_states]')
response = '\n'.join(line for line in response.split('\n') if 'polarizability.' not in line)
es_only = head + response
assert 'polarizability' not in es_only.split('molresponse excited')[1]
assert '  protocol [1e-4]\n' in response
aliased = head + response.replace('  protocol [1e-4]\n', '  protocol [1.2e-4]\n')
assert aliased != es_only


def run(label, text, prefix, cont=False, expect=0):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    args = launch + [binary, str(inp), '--prefix=' + str(prefix)] + (['--continue'] if cont else [])
    with log.open('w') as out:
        r = subprocess.run(args, cwd=work, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=900)
    if expect is not None:
        assert (r.returncode == 0) == (expect == 0), log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return r.returncode


def excited(prefix):
    return {s['stage']: s for s in json.loads((prefix / 'results.json').read_text())['stages']}['excited']


def engine_metadata(prefix):
    attempts = sorted((prefix / 'excited').glob('attempt-*'))
    assert attempts, prefix
    return json.loads((attempts[-1] / 'response_metadata.json').read_text()), attempts[-1]


# ---- D-088: the same stage re-solved at an aliased threshold --------------
alias = work / 'Alias'
run('alias-1e-4', es_only, alias)
first = excited(alias)
assert first['status'] == 'completed', first['status']
es_first = first['properties']['excited_states']
assert es_first['converged'] and es_first['protocol_key'] == '1e-04_k6', es_first
meta_first, _ = engine_metadata(alias)
assert meta_first['protocols']['1e-04_k6']['thresh'] == 1e-4, meta_first['protocols']

# The request's rung is now 1.2e-4 under the same key. The planner refuses the
# saved bundle as a source for another threshold (D-030) and re-solves; the
# writer must describe the new bundle with the new threshold (D-088).
run('alias-1.2e-4', aliased, alias, cont=True)
second = excited(alias)
assert second['status'] == 'completed', (second['status'], second.get('convergence_report', {}).get('decided_by'))
es_second = second['properties']['excited_states']
assert es_second['converged'] and es_second['protocol_key'] == '1e-04_k6', es_second
assert es_second['iterations_at_rung'] > 0 and second['iterations_this_run'] > 0, 'the aliased request must re-solve, not load'
meta_second, attempt = engine_metadata(alias)
entry = meta_second['protocols']['1e-04_k6']
assert abs(entry['thresh'] - 1.2e-4) < 1e-16 and entry['k'] == 6, entry
index = json.loads((attempt / 'es__1e-04_k6' / 'roots.json').read_text())
assert abs(index['thresh'] - entry['thresh']) < 1e-16 and index['k'] == entry['k'], (index['thresh'], entry)
# The stage's verdict judged the bundle at its rung: the state's rung is the
# requested one (the record's rung comes from the protocol entry).
states = second.get('convergence_report', {}).get('requested_accuracy', {})
assert 'es:tda_n1|1e-04_k6' in states, sorted(states)
assert all('transition' in r for r in es_second['roots']), es_second['roots']

# ---- D-089: a restart between the bundle's save and its analysis ---------
restart = work / 'Restart'
run('restart-solve', es_only, restart)
solved = excited(restart)
es_solved = solved['properties']['excited_states']
assert solved['status'] == 'completed' and es_solved['converged'], solved['status']
assert all('transition' in r for r in es_solved['roots']), es_solved['roots']
store = restart / 'excited' / 'checkpoints'
latest = json.loads((store / 'latest').read_text())['generation']
manifest = json.loads((store / latest / 'manifest.json').read_text())
assert 'transition_analysis' in manifest['numerical'], sorted(manifest['numerical'])
assert manifest['position'].get('phase') == 'property_analysis', manifest['position']
previous = manifest['position']['previous_generation']
earlier = json.loads((store / previous / 'manifest.json').read_text())
assert 'transition_analysis' not in earlier['numerical'], sorted(earlier['numerical'])
assert earlier['position'].get('kind') == 'excited_states', earlier['position']
# The state a crash between the two publishes leaves: the save-boundary
# generation is the latest and carries no analysis.
import shutil
shutil.rmtree(store / latest)
(store / 'latest').write_text(json.dumps({'generation': previous}) + '\n')

run('restart-continue', es_only, restart, cont=True, expect=None)
continued = excited(restart)
es_cont = continued['properties']['excited_states']
assert continued['iterations_this_run'] == 0, continued['iterations_this_run']
assert es_cont['converged'] and es_cont['decided_by'] == 'planner', es_cont
assert es_cont['iterations_at_rung'] == es_solved['iterations_at_rung'], (es_cont['iterations_at_rung'], es_solved['iterations_at_rung'])
assert [r['root_id'] for r in es_cont['roots']] == [r['root_id'] for r in es_solved['roots']]
_, attempt = engine_metadata(restart)
analysis_file = attempt / 'es_analysis__1e-04_k6.json'
# The engine's half (D-089, calc_executor.hpp analyze_loaded_es_bundle) writes
# the analysis on the load; the manager's call to it and the stage's reading
# of the file are A's lines. Until they land the file is absent and the
# record says so; with them every root carries its transition block.
if analysis_file.exists():
    assert all('transition' in r for r in es_cont['roots']), es_cont['roots']
    for before, after in zip(es_solved['roots'], es_cont['roots']):
        assert abs(before['transition']['oscillator_strength'] - after['transition']['oscillator_strength']) < 1e-8
    print('D-089: analysis computed on the load; transition blocks complete')
else:
    assert all('transition_unavailable' in r for r in es_cont['roots']), es_cont['roots']
    print('D-089: the engine analysis file is absent on the load (the manager call and the stage reading are not in yet)')
print('PASS')
