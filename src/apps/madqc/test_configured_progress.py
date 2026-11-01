#!/usr/bin/env python3
"""Explicit stall controls separate supplied-rung progress from accuracy."""
import argparse
import json
import sys
import os
from pathlib import Path
import re
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
work = named_workdir.work_directory('configured-progress-')
print('artifacts:', work, flush=True)
deck = (Path(args.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
# Experimental fixtures: exact-zero remains unmet; the supplied near-one
# improvement requirement makes the small study reach a stall promptly.
# progress.* acts on written targets, so the targets decide (DESIGN B8, Q5);
# the kinetic target is the only one, as dconv and econv are not written.
deck = deck.replace('  dconv 1e-5\n', '').replace('  econv 1e-6\n', '')
deck = deck.replace('moldft gs\n', 'moldft gs\n  convergence.use_default_criteria false\n'
                    '  accuracy.kinetic_trace_change 0\n'
                    '  progress.window 1\n  progress.minimum_relative_improvement 0.99999999\n')
deck = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter 12', deck)
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]
binary = str(Path(args.binary).resolve())

def run(label, text, stop_reason):
    inp = work / (label + '.in')
    inp.write_text(text)
    root = work / label
    with (work / (label + '.log')).open('w') as log:
        proc = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)],
                              cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=500)
    assert proc.returncode != 0  # unmet accuracy is never successful convergence
    stage = named_records.stages(root)[0]
    assert stage['status'] == 'unconverged'
    store = root / 'gs/checkpoints'
    manifest = json.loads((store / json.loads((store / 'latest').read_text())['generation'] / 'manifest.json').read_text())
    # B5: the stop reason names the exit (stalled or iteration_limit).
    assert manifest['stop_reason'] == stop_reason and manifest['archives'], manifest['stop_reason']
    report = stage['convergence_report']
    assert report['requested_accuracy']['status'] == 'unmet'
    return stage, report

stage, report = run('Stalled', deck, 'stalled')
assert stage['iterations_this_run'] == 6
history = report['iteration_history']
assert len({row['kinetic_trace']['sample']['threshold'] for row in history}) == 2
applied = [row['configured_progress']['action'] for row in history if row['configured_progress']['action_applied']]
assert applied == ['advance_supplied_rung', 'retain_unmet'], applied
assert report['configured_progress']['action'] == 'retain_unmet'

# An altered recorded decision is not accepted as valid continuation evidence.
invalid = work / 'InvalidProgress'
shutil.copytree(work / 'Stalled', invalid)
store = invalid / 'gs/checkpoints'
pointer = (store / 'latest').read_text()
manifest_file = store / json.loads(pointer)['generation'] / 'manifest.json'
manifest = json.loads(manifest_file.read_text())
report = manifest['metrics']['convergence_report']
for entry in report['iteration_history'] + [report]:
    has_next = entry['configured_progress']['has_next_supplied_rung']
    entry['configured_progress'] = {'status': 'disabled', 'action': 'continue',
                                    'has_next_supplied_rung': has_next, 'action_applied': False}
manifest_file.write_text(json.dumps(manifest))
with (work / 'invalid-progress.log').open('w') as log:
    proc = subprocess.run(launch + [binary, str(work / 'Stalled.in'), '--prefix=' + str(invalid), '--continue'],
                          cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=60)
assert proc.returncode != 0 and (store / 'latest').read_text() == pointer
assert 'invalid ground-state measurement' in (work / 'invalid-progress.log').read_text()

stage, report = run('Forced', deck.replace('moldft gs\n', 'moldft gs\n  diagnostic.iterations_per_rung 5\n'),
                    'iteration_limit')
assert stage['iterations_this_run'] == 10
assert not any(row['configured_progress']['action_applied'] for row in report['iteration_history'])
print('configured progress acceptance passed', flush=True)
