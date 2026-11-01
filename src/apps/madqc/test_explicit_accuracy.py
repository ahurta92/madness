#!/usr/bin/env python3
"""Explicit SCF criteria, saved-state reuse, unmet continuation and corruption."""
import argparse
import json
import sys
import os
import re
from pathlib import Path
import shutil
import subprocess
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)
import named_workdir  # noqa: E402  (D-108: removed when the test passes)

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--source', required=True)
parser.add_argument('--ranks', type=int, default=1)
parser.add_argument('--mpiexec', default='mpiexec')
args = parser.parse_args()
work = named_workdir.work_directory('explicit-accuracy-')
print('artifacts:', work, flush=True)
binary = str(Path(args.binary).resolve())
deck = (Path(args.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
# Deliberately loose experimental targets test opt-in plumbing, not
# calibration. dconv and econv are the density- and energy-change targets
# (DESIGN B8); here the written targets decide.
deck = deck.replace('dconv 1e-5', 'dconv 1000').replace('econv 1e-6', 'econv 1000')
controls = '  convergence.use_default_criteria false\n' + ''.join(
    '  accuracy.' + name + ' 1000\n' for name in ('bsh_residual', 'kinetic_trace_change'))
deck = deck.replace('moldft gs\n', 'moldft gs\n' + controls)
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
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
    return named_records.stages(root)[0]

def latest(root):
    store = root / 'gs/checkpoints'
    pointer = (store / 'latest').read_text()
    manifest = store / json.loads(pointer)['generation'] / 'manifest.json'
    return store, pointer, manifest

fresh = work / 'Fresh'
result = run('fresh', fresh)
report = result['convergence_report']
assert report['decided_by'] == 'requested_criteria'
assert report['requested_accuracy']['status'] == 'satisfied'
assert len(report['requested_accuracy']['criteria']) == 4
assert 'charge_density_change' in report['requested_accuracy']['criteria']
assert 'density_change' not in report['requested_accuracy']['criteria']
assert report['saved_state_validation']['explicit_accuracy']['status'] == 'satisfied'
repeat = run('repeat', fresh, resume=True)
assert repeat['iterations_this_run'] == 0 and repeat['state_id'] == result['state_id']

# Exact reuse must validate new configuration even when it would solve zero iterations.
for label, extra, base in (
        ('negative-window', '  progress.window -1\n', deck),
        ('missing-improvement', '  progress.window 1\n', deck),
        ('missing-targets', '  progress.window 1\n  progress.minimum_relative_improvement 0.1\n', deck.replace(controls, ''))):
    root = work / ('Invalid-' + label)
    shutil.copytree(fresh, root)
    store, pointer, _ = latest(root)
    run(label, root, text=base.replace('moldft gs\n', 'moldft gs\n' + extra), resume=True, success=False)
    expected = ('progress.window must be nonnegative' if label == 'negative-window'
                else 'progress.* stall rules act on your written targets' if label == 'missing-targets'
                else 'progress requires explicit accuracy targets')
    assert expected in (work / (label + '.log')).read_text()
    assert (store / 'latest').read_text() == pointer

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
        report['saved_state_validation']['explicit_accuracy']['criteria']['absolute_kinetic_trace_change']['status'] = 'unmet'
    manifest_file.write_text(json.dumps(manifest))
    run('invalid-' + location, root, resume=True, success=False)
    assert 'invalid ground-state measurement' in (work / ('invalid-' + location + '.log')).read_text()
    assert (store / 'latest').read_text() == pointer

# When the engine check decides, loose written targets are only reported:
# they cannot turn an unfinished engine into convergence (B8).
root = work / 'EngineUnmet'
short = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter 1',
               deck.replace('  convergence.use_default_criteria false\n', ''))
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
# maxiter is per run (B5, D-007): three more iterations.
tight = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter 3', tight)
result = run('tightened', root, text=tight, resume=True, success=False)
assert result['work'] == 'resumed' and result['iterations_this_run'] > 0
assert result['status'] == 'unconverged'
assert result['convergence_report']['requested_accuracy']['status'] == 'unmet'
_, pointer, manifest_file = latest(root)
assert pointer != original_pointer
manifest = json.loads(manifest_file.read_text())
assert manifest['stop_reason'] == 'iteration_limit'  # B5: the cap ended it
assert manifest['archives']  # usable orbitals retained in an immutable generation
print('explicit accuracy acceptance passed', flush=True)
