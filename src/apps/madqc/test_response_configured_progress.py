#!/usr/bin/env python3
"""FD configured stall decisions use only supplied rungs and explicit targets."""
import argparse
import json
import sys
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--source', required=True)
parser.add_argument('--ranks', type=int, default=1)
parser.add_argument('--mpiexec', default='mpiexec')
args = parser.parse_args()
work = Path(tempfile.mkdtemp(prefix='response-configured-progress-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
binary = str(Path(args.binary).resolve())
deck = (Path(args.source) / 'named_h2.in').read_text().replace(
    'polarizability.frequencies [0.0, 0.02]', 'polarizability.frequencies [0.0]')
# Deliberate stall fixture; these values are not calibrated scientific defaults.
# progress.* acts on written targets, so the targets decide (DESIGN B8); the
# kinetic target is the only one, as the response dconv is not written.
assert deck.count('  dconv 1e-4\n') == 1  # the response stage's dconv
deck = deck.replace('  dconv 1e-4\n', '')
deck = deck.replace('molresponse optical\n', 'molresponse optical\n'
                    '  convergence.use_default_criteria false\n'
                    '  accuracy.kinetic_trace_change 0\n  progress.window 1\n'
                    '  progress.minimum_relative_improvement 0.99999999\n')
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]

def run(label, root, text, resume=False):
    inp = work / (label + '.in')
    inp.write_text(text)
    with (work / (label + '.log')).open('w') as log:
        proc = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                              (['--continue'] if resume else []), cwd=work, env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=600)
    assert proc.returncode != 0
    return proc

def latest(root):
    store = root / 'optical/checkpoints'
    pointer = (store / 'latest').read_text()
    path = store / json.loads(pointer)['generation'] / 'manifest.json'
    return store, pointer, path, json.loads(path.read_text())

def entries(manifest):
    for channel in manifest['numerical']['engine_metadata']['fd_states'].values():
        for rung in channel.values():
            yield from rung.values()

root = work / 'Stalled'
run('stalled', root, deck)
result = named_records.stages(root)[-1]
assert result['status'] == 'unconverged' and result['iterations_this_run'] == 6
_, _, _, manifest = latest(root)
assert manifest['archives'] and not manifest['properties']
rows = sorted(entries(manifest), key=lambda e: -e['function_bsh_history'][0]['sample']['threshold'])
assert [e['configured_progress']['action'] for e in rows] == ['advance_supplied_rung', 'retain_unmet']
assert all(e['configured_progress']['action_applied'] and e['stalled'] and
           not e['converged'] and not e['accepted'] for e in rows)

for corruption in ('decision', 'fractional-window', 'missing-rung', 'inserted-rung', 'missing-progress', 'missing-registry'):
    bad = work / ('BadEvidence-' + corruption)
    shutil.copytree(root, bad)
    store, pointer, path, manifest = latest(bad)
    for entry in entries(manifest):
        if corruption == 'decision':
            for decision in entry['configured_progress_history'] + [entry['configured_progress']]:
                decision['status'] = 'disabled'
                decision['action'] = 'continue'
                decision['action_applied'] = False
        elif corruption == 'fractional-window':
            entry['checkpoint_progress']['window'] = 1.5
        elif corruption == 'missing-rung':
            entry['checkpoint_progress']['protocol'] = [1e-8]
        elif corruption == 'inserted-rung':
            ladder = entry['checkpoint_progress']['protocol']
            ladder.insert(1, (ladder[0] * ladder[1]) ** 0.5)
        else:
            for key in ('checkpoint_progress', 'configured_progress', 'configured_progress_history'):
                del entry[key]
    if corruption == 'missing-registry':
        del manifest['numerical']['engine_metadata']['checkpoint_requests']
    path.write_text(json.dumps(manifest))
    run('bad-' + corruption, bad, deck, resume=True)
    assert (store / 'latest').read_text() == pointer
    assert 'invalid response state' in (work / ('bad-' + corruption + '.log')).read_text()

# The final rung stalled with the request unmet, so --continue loads it and
# reports it unmet instead of stalling again (DESIGN B6 R2, D-005.1). The new
# generation keeps the originating request beside the current attempt's.
continued = work / 'Continued'
shutil.copytree(root, continued)
run('continued', continued, deck, resume=True)
result = named_records.stages(continued)[-1]
assert result['status'] == 'unconverged' and result['iterations_this_run'] == 0
states = result['continuation']['states']
assert [s['action'] for s in states] == ['loaded'] and states[0]['unmet']
_, _, _, manifest = latest(continued)
assert len(manifest['numerical']['engine_metadata']['checkpoint_requests']) >= 2
assert not manifest['properties']

forced = work / 'Forced'
forced_deck = deck.replace('molresponse optical\n', 'molresponse optical\n  diagnostic.iterations_per_rung 5\n')
run('forced', forced, forced_deck)
result = named_records.stages(forced)[-1]
assert result['status'] == 'unconverged' and result['iterations_this_run'] == 10
_, _, _, manifest = latest(forced)
assert not manifest['properties']
for entry in entries(manifest):
    assert entry['rung_exhausted'] and not entry['stalled']
    assert not any(d['action_applied'] for d in entry['configured_progress_history'])
    assert entry['checkpoint_progress']['minimum_relative_improvement'] == 0.99999999

store, pointer, _, _ = latest(root)
for label, text in [('zero-fraction', deck.replace('0.99999999', '0')),
                    ('one-fraction', deck.replace('0.99999999', '1')),
                    ('negative-window', deck.replace('progress.window 1', 'progress.window -1')),
                    ('no-target', deck.replace('  accuracy.kinetic_trace_change 0\n', ''))]:
    run(label, root, text, resume=True)
    assert (store / 'latest').read_text() == pointer
    # Without a target the targets cannot decide (B8): rejected at parse.
    expected = 'needs at least one written target' if label == 'no-target' else 'progress'
    assert expected in (work / (label + '.log')).read_text()
print('response configured progress acceptance passed', flush=True)
