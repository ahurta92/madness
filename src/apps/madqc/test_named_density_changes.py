#!/usr/bin/env python3
"""Report charge/spin changes without inventing new convergence thresholds."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='named-density-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
# This verified deck explicitly sets protocol, dconv, econv and maxiter.
# The new measurements have no stopping targets; no defaults are being tested.
deck = (Path(a.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
env = dict(os.environ, OPENBLAS_NUM_THREADS='1',
           OMP_NUM_THREADS='1', OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]

def run(label, root, text=deck, resume=False, success=True):
    inp = work / (label + '.in'); inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=240)
    assert (r.returncode == 0) == success, log.read_text()[-4000:]
    result = (root / 'results.json').read_text()
    (work / (label + '-results.json')).write_text(result)
    print(label, 'exit', r.returncode, flush=True)
    return json.loads(result)['stages'][0]

def check(measurements):
    charge = measurements['charge_density_change']
    spin = measurements['spin_density_change']
    for m in (charge, spin):
        assert m['active_gate'] is False
        assert m['target'] is None and m['effective_limit'] is None and m['satisfied'] is None
        assert m['units'] == 'a0^(-3/2)'
        assert m['reduction'] == 'L2_norm_of_field_difference'
        assert m['availability'] == measurements['density_change']['availability']
    assert charge['density_definition'] == 'rho_alpha + rho_beta'
    assert spin['density_definition'] == 'rho_alpha - rho_beta'
    if charge['availability'] == 'measured':
        assert math.isclose(charge['value'], 2 * measurements['density_change']['value'],
                            rel_tol=1e-9, abs_tol=1e-12)
        assert spin['value'] < 1e-12
    else:
        assert charge['availability'] == 'unavailable'
        assert charge['value'] is None and spin['value'] is None

root = work / 'Fresh'
fresh = run('fresh', root)
assert fresh['iterations_this_run'] == 8
report = fresh['convergence_report']
check(report['measurements'])
check(report['saved_state_validation']['measurements'])
for mf in (root / 'gs/checkpoints').glob('*/manifest.json'):
    m = json.loads(mf.read_text())
    if 'iteration' in m['position']:
        check(m['metrics']['convergence_report']['measurements'])
        if m['position']['iteration'] == 0:
            assert m['metrics']['convergence_report']['measurements']['charge_density_change']['value'] is None
repeat = run('reuse', root, resume=True)
assert repeat['iterations_this_run'] == 0 and repeat['convergence_report'] == report
partial = work / 'Partial'
short = deck.replace('maxiter 30', 'maxiter 1')
stopped = run('unfinished', partial, short, success=False)
check(stopped['convergence_report']['measurements'])
check(stopped['convergence_report']['saved_state_validation']['measurements'])
resumed = run('resume', partial, resume=True)
assert resumed['status'] == 'completed' and resumed['iterations_this_run'] > 0
check(resumed['convergence_report']['measurements'])
check(resumed['convergence_report']['saved_state_validation']['measurements'])
cases = [
    ('historical-charge', False, 'charge_density_change', {'value': -1}),
    ('saved-spin', True, 'spin_density_change', {'value': -1}),
    ('forbidden-target', True, 'charge_density_change', {'target': 1e-5}),
    ('forbidden-gate', False, 'spin_density_change', {'active_gate': True}),
]
for label, saved, name, edits in cases:
    invalid = work / ('Invalid-' + label); shutil.copytree(root, invalid)
    store = invalid / 'gs/checkpoints'; latest = (store / 'latest').read_text()
    mf = store / json.loads(latest)['generation'] / 'manifest.json'
    m = json.loads(mf.read_text())
    r = m['metrics']['convergence_report']
    if saved:
        r = r['saved_state_validation']
    r['measurements'][name].update(edits)
    mf.write_text(json.dumps(m))
    run('invalid-' + label, invalid, resume=True, success=False)
    assert 'invalid ground-state measurement' in (work / ('invalid-' + label + '.log')).read_text()
    assert (store / 'latest').read_text() == latest
# Optional old diagnostics stay absent on exact reuse; no invented zeros.
legacy = work / 'Legacy'; shutil.copytree(root, legacy)
store = legacy / 'gs/checkpoints'
mf = store / json.loads((store / 'latest').read_text())['generation'] / 'manifest.json'
m = json.loads(mf.read_text())
# Written before S4a too: no state record (numerical.states).
m['numerical'].pop('states')
r = m['metrics']['convergence_report']
for measurements in [r['measurements'], r['saved_state_validation']['measurements']]:
    del measurements['charge_density_change'], measurements['spin_density_change']
# Generations from before the charge-density measurement evaluated no written
# targets either (a written dconv targets that measurement now, DESIGN B8):
# their recorded evaluation is the empty one.
not_requested = {'status': 'not_requested', 'criteria': {}, 'comparison_relation': 'less_than_or_equal'}
for entry in [r, r['saved_state_validation']] + r['iteration_history']:
    entry['explicit_accuracy'] = dict(not_requested)
mf.write_text(json.dumps(m))
old = run('legacy', legacy, resume=True)
assert old['iterations_this_run'] == 0
assert 'charge_density_change' not in old['convergence_report']['measurements']
assert 'spin_density_change' not in old['convergence_report']['saved_state_validation']['measurements']
print('SCF density reporting passed', flush=True)
