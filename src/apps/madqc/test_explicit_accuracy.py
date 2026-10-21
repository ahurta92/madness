#!/usr/bin/env python3
"""Explicit SCF criteria, saved-state reuse, unmet continuation and corruption."""
import argparse
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--source', required=True)
parser.add_argument('--ranks', type=int, default=1)
parser.add_argument('--mpiexec', default='mpiexec')
args = parser.parse_args()
work = Path(tempfile.mkdtemp(prefix='explicit-accuracy-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
binary = str(Path(args.binary).resolve())
deck = (Path(args.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
# Deliberately loose experimental controls test opt-in plumbing, not calibration.
controls = ''.join('  accuracy.' + name + ' 1000\n' for name in
                   ('density_change', 'bsh_residual', 'energy_change', 'kinetic_trace_change'))
deck = deck.replace('moldft gs\n', 'moldft gs\n' + controls)
env = dict(os.environ, MAD_NUM_THREADS='2', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]

def run(label, root, text=deck, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    with (work / (label + '.log')).open('w') as output:
        result = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                                (['--continue'] if resume else []), cwd=work, env=env,
                                stdout=output, stderr=subprocess.STDOUT, timeout=600)
    assert (result.returncode == 0) == success, (work / (label + '.log')).read_text()[-5000:]
    print(label, 'exit', result.returncode, flush=True)
    return json.loads((root / 'results.json').read_text())['stages'][0]

def latest(root):
    store = root / 'gs/checkpoints'
    pointer = (store / 'latest').read_text()
    manifest = store / json.loads(pointer)['generation'] / 'manifest.json'
    return store, pointer, manifest

fresh = work / 'Fresh'
result = run('fresh', fresh)
report = result['convergence_report']
assert report['gate_policy'] == 'existing_engine_and_explicit_accuracy'
assert report['requested_accuracy']['status'] == 'satisfied'
assert len(report['requested_accuracy']['criteria']) == 4
assert 'charge_density_change' in report['requested_accuracy']['criteria']
assert 'density_change' not in report['requested_accuracy']['criteria']
assert report['saved_state_validation']['explicit_accuracy']['status'] == 'satisfied'
repeat = run('repeat', fresh, resume=True)
assert repeat['iterations_this_run'] == 0 and repeat['state_id'] == result['state_id']

for location in ('input_phase', 'history_phase', 'saved_phase', 'accuracy_flag',
                 'comparison_threshold', 'comparison_k', 'comparison_run_id',
                 'comparison_stage', 'comparison_iteration', 'comparison_trace',
                 'history_stage', 'history_iteration', 'history_density', 'history_bsh'):
    root = work / ('Invalid-' + location)
    shutil.copytree(fresh, root)
    store, pointer, manifest_file = latest(root)
    manifest = json.loads(manifest_file.read_text())
    report = manifest['metrics']['convergence_report']
    if location == 'input_phase':
        manifest['metrics']['kinetic_trace']['phase'] = 'incompatible_phase'
    elif location == 'history_phase':
        report['iteration_history'][0]['kinetic_trace']['phase'] = 'incompatible_phase'
    elif location == 'saved_phase':
        report['saved_state_validation']['kinetic_trace']['phase'] = 'incompatible_phase'
    elif location.startswith('comparison_'):
        comparison = report['saved_state_validation']['comparison']
        field = location.removeprefix('comparison_')
        if field == 'trace':
            comparison['kinetic_trace']['value'] += 1
        elif field in ('stage', 'run_id'):
            comparison[field] += '-wrong'
        else:
            comparison[field] += 1
    elif location.startswith('history_'):
        entry = report['iteration_history'][1]
        field = location.removeprefix('history_')
        if field == 'stage':
            entry['kinetic_trace']['sample']['stage'] += '-wrong'
        elif field == 'iteration':
            entry['kinetic_trace']['sample']['iteration'] += 2
        elif field == 'density':
            entry['measurements']['spin_density_change']['availability'] = 'invalid'
        else:
            entry['measurements']['bsh_residual']['per_function']['blocks']['alpha'][0]['availability'] = 'invalid'
    else:
        report['saved_state_validation']['explicit_accuracy']['criteria']['delta_kinetic_energy']['status'] = 'unmet'
    manifest_file.write_text(json.dumps(manifest))
    run('invalid-' + location, root, resume=True, success=False)
    assert 'invalid ground-state measurement' in (work / ('invalid-' + location + '.log')).read_text()
    assert (store / 'latest').read_text() == pointer

# Loose explicit criteria cannot turn an unfinished engine into convergence.
root = work / 'EngineUnmet'
short = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter 1', deck)
result = run('engine-unmet', root, text=short, success=False)
assert result['status'] == 'unconverged'
assert result['convergence_report']['requested_accuracy']['status'] == 'satisfied'
_, _, manifest_file = latest(root)
assert json.loads(manifest_file.read_text())['stop_reason'] == 'iteration_limit'

root = work / 'Tightened'
shutil.copytree(fresh, root)
_, original_pointer, _ = latest(root)
tight = deck.replace('accuracy.kinetic_trace_change 1000', 'accuracy.kinetic_trace_change 0')
# No universal stopping tolerance is inferred: exact-zero is deliberately unmet.
tight = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter 3', tight)
result = run('tightened', root, text=tight, resume=True, success=False)
assert result['work'] == 'resumed' and result['iterations_this_run'] > 0
assert result['status'] == 'unconverged'
assert result['convergence_report']['requested_accuracy']['status'] == 'unmet'
_, pointer, manifest_file = latest(root)
assert pointer != original_pointer
manifest = json.loads(manifest_file.read_text())
assert manifest['stop_reason'] == 'unmet_requested_accuracy'
assert manifest['archives']  # usable orbitals retained in an immutable generation
print('explicit accuracy acceptance passed', flush=True)
