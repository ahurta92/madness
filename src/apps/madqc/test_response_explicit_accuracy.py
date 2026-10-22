#!/usr/bin/env python3
"""Named FD criteria must certify saved outputs before property eligibility."""
import argparse
import json
import os
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
work = Path(tempfile.mkdtemp(prefix='response-explicit-accuracy-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
binary = str(Path(args.binary).resolve())
deck = (Path(args.source) / 'named_h2.in').read_text()
deck = deck.replace('polarizability.frequencies [0.0, 0.02]', 'polarizability.frequencies [0.0]')
# Loose values are acceptance fixtures, never suggested scientific tolerances.
# The written targets decide (DESIGN B8); the response dconv is the
# induced-density target.
assert deck.count('  dconv 1e-4\n') == 1  # the response stage's dconv
deck = deck.replace('  dconv 1e-4\n', '  dconv 1000\n')
controls = '  convergence.use_default_criteria false\n' + ''.join(
    '  accuracy.' + metric + ' 1000\n' for metric in ('bsh_residual', 'kinetic_trace_change'))
assert deck.count('molresponse optical\n') == 1
deck = deck.replace('molresponse optical\n', 'molresponse optical\n' + controls)
assert 'polarizability.frequencies [0.0]' in deck and controls in deck
env = dict(os.environ, MAD_NUM_THREADS='2', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]

def run(label, root, text, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as output:
        proc = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                              (['--continue'] if resume else []), cwd=work, env=env,
                              stdout=output, stderr=subprocess.STDOUT, timeout=600)
    assert (proc.returncode == 0) == success, log.read_text()[-5000:]
    print(label, 'exit', proc.returncode, flush=True)
    return json.loads((root / 'results.json').read_text())['stages']

def latest(root):
    store = root / 'optical/checkpoints'
    pointer = (store / 'latest').read_text()
    manifest_path = store / json.loads(pointer)['generation'] / 'manifest.json'
    return store, pointer, manifest_path, json.loads(manifest_path.read_text())

def entries(manifest):
    for channel in manifest['numerical']['engine_metadata']['fd_states'].values():
        for rung in channel.values():
            yield from rung.values()

fresh = work / 'Fresh'
stages = run('fresh', fresh, deck)
assert stages[-1]['status'] == 'completed'
_, pointer, _, manifest = latest(fresh)
assert len(manifest['numerical']['request']['accuracy_targets']) == 3
checked = 0
for entry in entries(manifest):
    if not entry['converged']:
        continue
    checked += 1
    assert not entry['accepted']
    assert entry['explicit_accuracy']['status'] == 'satisfied'
    assert entry['saved_state_validation']['explicit_accuracy']['status'] == 'satisfied'
    criterion = entry['saved_state_validation']['explicit_accuracy']['criteria']['bsh_residual']
    raw = entry['saved_state_validation']['measurements']['bsh_residual']['per_function']['blocks']
    assert criterion['value'] == max(f['value'] for block in raw.values() for f in block)
    assert criterion['target']['quantity'] == 'maximum_raw_fixed_point_function_L2'
assert checked > 0
repeat = run('repeat', fresh, deck, resume=True)
assert repeat[-1]['iterations_this_run'] == 0
assert latest(fresh)[1] == pointer

# Derived scientific decisions and their comparison identities must be checked.
for corruption in ('decision', 'saved_value', 'history_phase', 'history_index', 'saved_comparison', 'input_density_phase', 'input_density_sample',
                   'input_bsh_sample', 'missing_targets'):
    root = work / ('Corrupt-' + corruption)
    shutil.copytree(fresh, root)
    store, before, path, damaged = latest(root)
    entry = next(entries(damaged))
    if corruption == 'decision':
        entry['explicit_accuracy']['status'] = 'unmet'
    elif corruption == 'saved_value':
        entry['saved_state_validation']['explicit_accuracy']['criteria']['bsh_residual']['value'] += 1
    elif corruption == 'history_phase':
        entry['function_bsh_history'][0]['phase'] = 'saved_response_before_archive'
    elif corruption == 'history_index':
        entry['function_kinetic_history'][0]['sample']['iteration'] = 99
    elif corruption == 'saved_comparison':
        entry['saved_state_validation']['comparison']['threshold'] *= 2
    elif corruption == 'input_density_phase':
        entry['convergence_measurements']['density_change']['phase'] = 'saved_response_before_archive'
    elif corruption == 'input_density_sample':
        entry['convergence_measurements']['density_change']['sample']['iteration'] += 1
    elif corruption == 'input_bsh_sample':
        entry['convergence_measurements']['bsh_residual']['sample']['attempt_directory'] += '-other'
    else:
        del entry['checkpoint_accuracy_targets']
    path.write_text(json.dumps(damaged))
    run('corrupt-' + corruption, root, deck, resume=True, success=False)
    assert 'invalid response state' in (work / ('corrupt-' + corruption + '.log')).read_text()
    assert (store / 'latest').read_text() == before

# A tighter request keeps compatible amplitudes but cannot inherit properties.
tight = work / 'Tight'
shutil.copytree(fresh, tight)
tight_deck = deck.replace('accuracy.kinetic_trace_change 1000', 'accuracy.kinetic_trace_change 0')
head, response = tight_deck.split('molresponse optical\n', 1)
tight_deck = head + 'molresponse optical\n' + response.replace('maxiter 30', 'maxiter 3')
assert 'accuracy.kinetic_trace_change 0' in tight_deck and 'maxiter 3' in tight_deck.split('molresponse optical\n')[1]
stages = run('tight', tight, tight_deck, resume=True, success=False)
assert stages[-1]['iterations_this_run'] > 0 and latest(tight)[1] != pointer
assert stages[-1]['status'] == 'unconverged'
_, _, _, manifest = latest(tight)
assert manifest['archives'] and manifest['stop_reason'] != 'converged'
assert not manifest['properties']
assert any(e['saved_state_validation']['explicit_accuracy']['status'] == 'unmet'
           for e in entries(manifest))

# Explicit diagnostic exhaustion advances rungs without accepted convergence.
forced = work / 'Forced'
shutil.copytree(fresh, forced)
forced_deck = deck.replace('accuracy.kinetic_trace_change 1000', 'accuracy.kinetic_trace_change 0').replace('molresponse optical\n', 'molresponse optical\n  diagnostic.iterations_per_rung 3\n')
forced_result = run('forced', forced, forced_deck, resume=True, success=False)
assert forced_result[-1]['iterations_this_run'] == 6
_, _, _, manifest = latest(forced)
assert not manifest['properties']
for entry in entries(manifest):
    assert entry['rung_exhausted'] and entry['iter'] == 3
    assert not entry['converged'] and not entry['accepted']
forced_again = run('forced-again', forced, forced_deck, resume=True, success=False)
assert forced_again[-1]['iterations_this_run'] == 6

# No local driven-response energy is fabricated from the source overlap.
unsupported = deck.replace('molresponse optical\n', 'molresponse optical\n  accuracy.energy_change 1\n')
run('unsupported-energy', fresh, unsupported, resume=True, success=False)
assert 'accuracy.energy_change is not applicable' in (work / 'unsupported-energy.log').read_text()
assert latest(fresh)[1] == pointer
# A driven Full response must retain both X and Y in its raw criteria.
full = work / 'Full'
full_deck = deck.replace('polarizability.frequencies [0.0]', 'polarizability.frequencies [0.02]')
full_result = run('full', full, full_deck)
assert full_result[-1]['status'] == 'completed'
_, _, _, full_manifest = latest(full)
for entry in entries(full_manifest):
    blocks = entry['saved_state_validation']['measurements']['bsh_residual']['per_function']['blocks']
    assert 'x_alpha' in blocks and 'y_alpha' in blocks
    if entry['converged']:
        assert entry['saved_state_validation']['explicit_accuracy']['status'] == 'satisfied'
print('response explicit accuracy acceptance passed', flush=True)
