#!/usr/bin/env python3
"""Check scientific evidence survives publication and continuation, not just results."""
import argparse
import json
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
work = Path(tempfile.mkdtemp(prefix='named-convergence-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
binary = str(Path(a.binary).resolve())
deck = (Path(a.source) / 'named_h2.in').read_text()
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]

def run(label, root, text=deck, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=360)
    assert (r.returncode == 0) == success, log.read_text()[-5000:]
    print(label, 'exit', r.returncode, flush=True)
    return json.loads((root / 'results.json').read_text())['stages']

def manifests(root, stage):
    return [json.loads(f.read_text()) for f in (root / stage / 'checkpoints').glob('*/manifest.json')]

fresh = work / 'Fresh'
stages = run('fresh', fresh)
gs, response = stages
assert 'convergence_report' in gs, 'ground convergence measurement report missing'
g = gs['convergence_report']
assert set(g['measurements']) == {'density_change', 'bsh_residual', 'energy_change', 'delta_kinetic_energy'}
validation = g['saved_state_validation']
assert validation['state_id'] == gs['state_id']
assert validation['affects_stopping'] is False
for name, m in validation['measurements'].items():
    assert m['availability'] == 'measured', (name, m)
    assert m['valid_for_saved_state'] is True
    assert m['value'] >= 0
assert validation['comparison']['iteration'] < g['measurements']['energy_change']['sample']['iteration']
assert g['gate_policy'] == 'existing_engine'
assert g['state_id'] == gs['state_id']
for name, m in g['measurements'].items():
    assert m['availability'] == 'measured', (name, m)
    assert m['value'] >= 0 and m['phase']
    assert m['sample']['iteration'] >= 0
    assert m['valid_for_saved_state'] is False
assert g['measurements']['energy_change']['active_gate'] is False
assert g['measurements']['energy_change']['target'] == 1e-6
assert g['measurements']['delta_kinetic_energy']['active_gate'] is False
r = response['convergence_report']
assert r['reference']['state_id'] == gs['state_id']
assert r['reference']['energy_change'] == g['measurements']['energy_change']
assert r['reference']['saved_state_validation'] == g['saved_state_validation']
assert len(r['channels']) == 4  # static + dynamic, two supplied rungs
for c in r['channels']:
    assert c['state_id'] == response['state_id'] and c['archive']
    m = c['measurements']
    assert m['density_change']['availability'] == 'measured'
    assert m['bsh_residual']['availability'] == 'measured'
    assert m['density_change']['effective_limit'] > 0
    assert m['bsh_residual']['valid_for_saved_state'] is False
    assert m['energy_change']['availability'] == 'not_applicable'
    assert m['delta_kinetic_energy']['availability'] == 'not_applicable'
    kinetics = c['function_kinetic_history']
    assert len(kinetics) == c['iterations_at_save'] + 1
    assert [k['sample']['iteration'] for k in kinetics] == list(range(c['iterations_at_save'] + 1))
    blocks = {'x_alpha'} if c['frequency'] == 0 else {'x_alpha', 'y_alpha'}
    for k in kinetics:
        assert k['reduction'] == 'none' and k['active_gate'] is False
        assert set(k['blocks']) == blocks
        for functions in k['blocks'].values():
            assert len(functions) == 1  # H2: one occupied spatial orbital
            assert functions[0]['function_index'] == 0
            assert functions[0]['availability'] == 'measured'
            assert functions[0]['value'] >= 0
    final = c['saved_state_validation']['measurements']['kinetic_energy']
    assert final['valid_for_saved_state'] is True
    for block in blocks:
        assert abs(final['blocks'][block][0]['value'] - kinetics[-1]['blocks'][block][0]['value']) < 1e-10
    assert c['outcome'] == 'converged'
    validation = c['saved_state_validation']
    assert validation['state_id'] == response['state_id']
    assert validation['affects_stopping'] is False
    for name in ('density_change', 'bsh_residual'):
        v = validation['measurements'][name]
        assert v['availability'] == 'measured' and v['valid_for_saved_state'] is True
        assert v['value'] >= 0 and v['effective_limit'] > 0
    assert validation['measurements']['delta_kinetic_energy']['availability'] == 'not_applicable'
# First differences at every SCF rung must be unavailable, never invented zeros.
first = [m for m in manifests(fresh, 'gs') if m['position'].get('iteration') == 0]
assert len(first) == 2
for m in first:
    evidence = m['metrics']['convergence_report']['measurements']
    for name in ('energy_change', 'density_change', 'delta_kinetic_energy'):
        assert evidence[name]['availability'] == 'unavailable' and evidence[name]['value'] is None
# Reported changes must agree with independently differenced saved samples.
by_rung = {}
for m in manifests(fresh, 'gs'):
    if 'iteration' in m['position']:
        by_rung.setdefault(m['numerical']['threshold'], []).append(m)
for samples in by_rung.values():
    samples.sort(key=lambda m: m['position']['iteration'])
    for before, after in zip(samples, samples[1:]):
        values = after['metrics']['convergence_report']['measurements']
        delta = abs(after['metrics']['energy']['value'] - before['metrics']['energy']['value'])
        assert abs(values['energy_change']['value'] - delta) < 1e-14
        delta_t = abs(after['metrics']['kinetic_energy']['value'] - before['metrics']['kinetic_energy']['value'])
        assert abs(values['delta_kinetic_energy']['value'] - delta_t) < 1e-14
# Final energy and kinetic differences must use the recorded prior iteration,
# never the almost-identical state just before final diagonalization.
v = g['saved_state_validation']
comparison = v['comparison']
previous = next(m for m in manifests(fresh, 'gs')
                if m['numerical']['threshold'] == comparison['threshold']
                and m['position'].get('iteration') == comparison['iteration'])
assert comparison['energy'] == previous['metrics']['energy']['value']
assert comparison['kinetic_energy'] == previous['metrics']['kinetic_energy']['value']
for metric, key in [('energy_change', 'energy'), ('delta_kinetic_energy', 'kinetic_energy')]:
    assert abs(v['measurements'][metric]['value'] - abs(v[key] - comparison[key])) < 1e-14
# Diagnostic evaluations must not add solver steps to the established H2 case.
assert gs['iterations_this_run'] == 8
assert response['iterations_this_run'] == 18
# Identical reuse must preserve recorded samples rather than fabricate fresh ones.
for stage in ('gs', 'optical'):
    (fresh / stage / 'results.json').unlink()
repeat = run('repeat', fresh, resume=True)
for before, after in zip(stages, repeat):
    assert after['iterations_this_run'] == 0
    assert before['convergence_report'] == after['convergence_report']
# One iteration creates a missing density comparison, not a zero change verdict.
partial = work / 'Partial'
shutil.copytree(fresh, partial)
shutil.rmtree(partial / 'optical')
short = 'maxiter 1'.join(deck.rsplit('maxiter 30', 1))
un = run('unfinished', partial, short, resume=True, success=False)[1]
assert un['status'] == 'unconverged'
assert all(c['saved_state_validation']['measurements']['density_change']['availability'] == 'measured'
           for c in un['convergence_report']['channels'])
assert all(c['outcome'] == 'unfinished' for c in un['convergence_report']['channels'])
assert all(c['measurements']['density_change']['value'] is None
           for c in un['convergence_report']['channels'])
resumed = run('resume', partial, resume=True)[1]
assert resumed['status'] == 'completed' and resumed['iterations_this_run'] > 0
# Independent check: an N=1 saved output is the input of step 2 in an N=2
# solve from the same initial reference. Its raw BSH and density differences
# must match the ordinary solver's next historical samples (static and full).
point_deck = 'protocol [1e-6]'.join(deck.rsplit('protocol [1e-4, 1e-6]', 1))
point_reports = []
for count in (1, 2):
    root = work / ('PointSteps' + str(count))
    shutil.copytree(fresh, root)
    shutil.rmtree(root / 'optical')
    text = ('maxiter ' + str(count)).join(point_deck.rsplit('maxiter 30', 1))
    point_reports.append(run('point-' + str(count), root, text, resume=True, success=False)[1])
for one, two in zip(point_reports[0]['convergence_report']['channels'],
                    point_reports[1]['convergence_report']['channels']):
    assert one['frequency'] == two['frequency']
    # Saved per-function values independently match the next actual iterate.
    saved = one['saved_state_validation']['measurements']['kinetic_energy']
    next_input = two['measurements']['kinetic_energy']
    for block in saved['blocks']:
        for a, b in zip(saved['blocks'][block], next_input['blocks'][block]):
            assert a['function_index'] == b['function_index']
            assert abs(a['value'] - b['value']) < 1e-10
    assert one['frequency'] == two['frequency']
    for metric in ('bsh_residual', 'density_change'):
        saved = one['saved_state_validation']['measurements'][metric]['value']
        next_input = two['measurements'][metric]['value']
        assert abs(saved - next_input) < 1e-10, (metric, saved, next_input)
# For SCF, step 3's historical BSH is step 2's input residual. A separate
# one-step solve archives precisely that input. Final rotation still compares
# to a genuine preceding iteration as checked above.
ground_deck = deck.split('molresponse optical')[0].replace('protocol [1e-4, 1e-6]', 'protocol [1e-6]')
ground_reports = []
for count in (1, 3):
    root = work / ('GroundSteps' + str(count))
    text = ground_deck.replace('maxiter 30', 'maxiter ' + str(count))
    ground_reports.append(run('ground-' + str(count), root, text, success=False)[0])
one = ground_reports[0]['convergence_report']['saved_state_validation']['measurements']
three = ground_reports[1]['convergence_report']['measurements']
assert abs(one['bsh_residual']['value'] - three['bsh_residual']['value']) < 1e-9
step2 = next(m for m in manifests(work / 'GroundSteps3', 'gs') if m['position'].get('iteration') == 1)
for metric in ('density_change', 'energy_change', 'delta_kinetic_energy'):
    expected = step2['metrics']['convergence_report']['measurements'][metric]['value']
    assert abs(one[metric]['value'] - expected) < 1e-9, (metric, one[metric]['value'], expected)
# A corrupt scientific verdict must not trigger an implicit new solve.
invalid = work / 'Invalid'
shutil.copytree(fresh, invalid)
store = invalid / 'optical/checkpoints'
latest = (store / 'latest').read_text()
mfile = store / json.loads(latest)['generation'] / 'manifest.json'
m = json.loads(mfile.read_text())
m['stop_reason'] = 'unfinished'
m['properties'] = []
meta = m['numerical']['engine_metadata']
meta['properties'] = {}
for channel in meta['fd_states'].values():
    for rung in channel.values():
        for e in rung.values():
            e.update(converged=False, diverged=True)
mfile.write_text(json.dumps(m))
run('invalid', invalid, resume=True, success=False)
assert 'invalid response state' in (work / 'invalid.log').read_text()
assert (store / 'latest').read_text() == latest
# Invalid ground evidence must not be reused or inherited by response either.
invalid_gs = work / 'InvalidGround'
shutil.copytree(fresh, invalid_gs)
store = invalid_gs / 'gs/checkpoints'
latest = (store / 'latest').read_text()
mfile = store / json.loads(latest)['generation'] / 'manifest.json'
m = json.loads(mfile.read_text())
m['metrics']['convergence_report']['measurements']['energy_change'].update(availability='invalid', value=None)
mfile.write_text(json.dumps(m))
run('invalid-ground', invalid_gs, resume=True, success=False)
assert 'invalid ground-state measurement' in (work / 'invalid-ground.log').read_text()
assert (store / 'latest').read_text() == latest
# Invalid revalidation evidence must fail before replacing either publication.
for stage in ('gs', 'optical'):
    root = work / ('InvalidValidation-' + stage)
    shutil.copytree(fresh, root)
    store = root / stage / 'checkpoints'
    latest = (store / 'latest').read_text()
    mfile = store / json.loads(latest)['generation'] / 'manifest.json'
    m = json.loads(mfile.read_text())
    if stage == 'gs':
        validation = m['metrics']['convergence_report']['saved_state_validation']
    else:
        channels = m['numerical']['engine_metadata']['fd_states']
        entry = next(iter(next(iter(next(iter(channels.values())).values())).values()))
        validation = entry['saved_state_validation']
    validation['measurements']['bsh_residual'].update(availability='invalid', value=None)
    mfile.write_text(json.dumps(m))
    run('invalid-validation-' + stage, root, resume=True, success=False)
    assert (store / 'latest').read_text() == latest
# A corrupt individual value must not hide behind a "measured" envelope.
for location in ('history', 'saved'):
    root = work / ('InvalidFunctionKinetic-' + location)
    shutil.copytree(fresh, root)
    store = root / 'optical/checkpoints'
    latest = (store / 'latest').read_text()
    mfile = store / json.loads(latest)['generation'] / 'manifest.json'
    m = json.loads(mfile.read_text())
    channels = m['numerical']['engine_metadata']['fd_states']
    entry = next(iter(next(iter(next(iter(channels.values())).values())).values()))
    record = (entry['function_kinetic_history'][0] if location == 'history' else
              entry['saved_state_validation']['measurements']['kinetic_energy'])
    record['blocks']['x_alpha'][0]['value'] = None
    mfile.write_text(json.dumps(m))
    run('invalid-function-kinetic-' + location, root, resume=True, success=False)
    assert 'invalid response state: per-function kinetic' in (work / ('invalid-function-kinetic-' + location + '.log')).read_text()
    assert (store / 'latest').read_text() == latest
# Old generations stay reusable but cannot acquire invented measurements.
old = work / 'Old'
shutil.copytree(fresh, old)
for stage in ('gs', 'optical'):
    store = old / stage / 'checkpoints'
    mfile = store / json.loads((store / 'latest').read_text())['generation'] / 'manifest.json'
    m = json.loads(mfile.read_text())
    if stage == 'gs':
        m['metrics'].pop('convergence_report')
    else:
        for channel in m['numerical']['engine_metadata']['fd_states'].values():
            for rung in channel.values():
                for e in rung.values():
                    e.pop('convergence_measurements')
                    e.pop('saved_state_validation', None)
                    e.pop('function_kinetic_history', None)
    mfile.write_text(json.dumps(m))
old_results = run('old', old, resume=True)
assert all(s['iterations_this_run'] == 0 for s in old_results)
assert old_results[0]['convergence_report']['measurements']['energy_change']['value'] is None
assert old_results[1]['convergence_report']['channels'][0]['measurements']['bsh_residual']['value'] is None
assert 'saved_state_validation' not in old_results[0]['convergence_report']
assert all('saved_state_validation' not in c and 'function_kinetic_history' not in c
           for c in old_results[1]['convergence_report']['channels'])
print('named convergence reporting, saved-state validation and continuation passed', flush=True)
