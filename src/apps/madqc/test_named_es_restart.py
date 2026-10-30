#!/usr/bin/env python3
"""Excited-state bundles across runs on one named stage (D-088, D-089).

D-088: a re-solve at an aliased threshold (1.2e-4 and 1e-4 share the key
1e-04_k6) redefines the protocol entry it is saved under, so the state record,
the engine's protocol entry and the bundle's own index describe the same rung
and the stage's verdict is on the requested rung. Before the fix the entry
kept the old threshold, same_rung failed the request and the stage ended
unconverged although the bundle had converged.
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
print('PASS')
