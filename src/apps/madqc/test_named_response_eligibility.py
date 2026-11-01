#!/usr/bin/env python3
"""A residual within tolerance cannot certify an unmeasured density change."""
import argparse
import json
import sys
import os
from pathlib import Path
import shutil
import subprocess
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)
import named_workdir  # noqa: E402  (D-108: removed when the test passes)

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = named_workdir.work_directory('response-eligibility-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1',
           OMP_NUM_THREADS='1', OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
deck = (Path(a.source) / 'named_h2.in').read_text()
ground, response = deck.split('molresponse optical', 1)
# Deliberately loose diagnostic target: both residual limits hold at step 1,
# but the first input has no preceding density. This is not an accuracy test.
response = response.replace('protocol [1e-4, 1e-6]', 'protocol [1e-6]')
response = response.replace('dconv 1e-4', 'dconv 10')

def run(label, root, text, success, resume=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        proc = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                              (['--continue'] if resume else []), cwd=work, env=env,
                              stdout=out, stderr=subprocess.STDOUT, timeout=240)
    assert (proc.returncode == 0) == success, (label, proc.returncode, log.read_text()[-4000:])
    print(label, 'exit', proc.returncode, flush=True)
    results = (root / 'results.json').read_text()
    (work / (label + '-results.json')).write_text(results)
    return named_records.stages(root)

def response_deck(count):
    return ground + 'molresponse optical' + response.replace('maxiter 30', 'maxiter ' + str(count))

def check_unfinished(stage):
    assert stage['status'] == 'unconverged', stage
    # Partial property estimates remain reportable, but cannot be certified.
    assert all(not value['converged']
               for value in stage['properties']['polarizability']['values'])
    channels = stage['convergence_report']['channels']
    assert {c['frequency'] for c in channels} == {0.0, 0.02}
    for c in channels:
        assert c['iterations_at_save'] == 1
        assert c['outcome'] == 'unfinished'
        assert c['measurements']['bsh_residual']['satisfied'] is True
        assert c['measurements']['density_change']['availability'] == 'unavailable'
        assert c['measurements']['density_change']['satisfied'] is None
    return channels

root = work / 'Run'
run('ground', root, ground, True, resume=False)
one = run('one-step', root, response_deck(1), False)[1]
check_unfinished(one)
# Restart supplies functions, not a second historical sample in this attempt.
one_again = run('one-step-restart', root, response_deck(1), False)[1]
check_unfinished(one_again)
assert one_again['iterations_this_run'] > 0
completed = run('two-step-restart', root, response_deck(2), True)[1]
assert completed['status'] == 'completed'
assert 'polarizability' in completed['properties']
channels = completed['convergence_report']['channels']
assert {c['frequency'] for c in channels} == {0.0, 0.02}
for c in channels:
    assert c['iterations_at_save'] == 2
    assert c['outcome'] == 'converged'
    assert c['measurements']['density_change']['availability'] == 'measured'
    assert c['measurements']['density_change']['satisfied'] is True
    assert c['measurements']['bsh_residual']['satisfied'] is True
reused = run('exact-reuse', root, response_deck(2), True)[1]
assert reused['iterations_this_run'] == 0
assert reused['convergence_report'] == completed['convergence_report']
# Old one-step certificates and stalled records must not bypass the corrected
# executor through completed checkpoint reuse. Only metadata is changed here;
# the actual archived functions and their integrity information stay intact.
for label, edits in [('invalid-one-step', {'iter': 1}),
                     ('invalid-stalled', {'stalled': True})]:
    invalid = work / label
    shutil.copytree(root, invalid)
    store = invalid / 'optical/checkpoints'
    latest = (store / 'latest').read_text()
    mfile = store / json.loads(latest)['generation'] / 'manifest.json'
    manifest = json.loads(mfile.read_text())
    states = manifest['numerical']['engine_metadata']['fd_states']
    entry = next(iter(next(iter(next(iter(states.values())).values())).values()))
    assert entry['converged'] and not entry['accepted']
    entry.update(edits)
    mfile.write_text(json.dumps(manifest))
    run(label, invalid, response_deck(2), False)
    assert 'ineligible strict convergence' in (work / (label + '.log')).read_text()
    assert (store / 'latest').read_text() == latest
print('response convergence eligibility passed', flush=True)
