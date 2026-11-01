#!/usr/bin/env python3
"""Check scientific evidence survives publication and continuation, not just
results. Runs at the default checkpoint cadence (DESIGN A3): only decision
points are published, so per-iteration evidence is read from the history the
stop manifest carries."""
import argparse
import json
import sys
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)

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
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
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
    return named_records.stages(root)

def strip_written_targets(entry):
    """Make a response entry look as one from before step 2c: the response
    dconv was not yet a written target (DESIGN B8), so such entries carry no
    written-target evidence."""
    for key in ('checkpoint_accuracy_targets', 'checkpoint_decided_by', 'explicit_accuracy',
                'explicit_accuracy_history', 'checkpoint_iteration_count',
                'checkpoint_diagnostic_iterations', 'rung_exhausted'):
        entry.pop(key, None)
    entry.get('saved_state_validation', {}).pop('explicit_accuracy', None)

def latest_manifest(root, stage):
    store = root / stage / 'checkpoints'
    return json.loads((store / json.loads((store / 'latest').read_text())['generation'] /
                       'manifest.json').read_text())


def gs_history(root, stage='gs'):
    """The per-iteration records the stop manifest carries (DESIGN A3: by
    default only decision points are published, and the history travels with
    them). Each row holds the measurements its sample's manifest would hold."""
    return latest_manifest(root, stage)['metrics']['convergence_report']['iteration_history']

def check_function_bsh(metric, reduction):
    record = metric['per_function']
    assert record['quantity'] == 'raw_fixed_point_function_L2'
    assert record['availability'] == 'measured'
    assert record['normalization'] == 'none' and record['reduction'] == 'none'
    assert record['active_gate'] is False
    values = []
    for functions in record['blocks'].values():
        for index, function in enumerate(functions):
            assert function['function_index'] == index
            assert function['availability'] == 'measured'
            assert math.isfinite(function['value']) and function['value'] >= 0
            values.append(function['value'])
    assert values
    expected = max(values) if reduction == 'max' else math.sqrt(sum(v*v for v in values))
    assert abs(metric['value'] - expected) < 1e-10

def check_bsh_comparison(metric):
    comparison = metric['maximum_comparison']
    values = [f['value'] for block in metric['per_function']['blocks'].values() for f in block]
    assert comparison['value'] == max(values)
    assert comparison['active_gate'] is False
    assert comparison['effective_limit'] == metric['effective_limit']
    assert comparison['satisfied'] == (max(values) <= metric['effective_limit'])
    assert math.isclose(metric['effective_limit'], metric['absolute_limit'] + metric['relative_factor'] * metric['normalization_scale'], rel_tol=1e-15, abs_tol=0.0)
    assert comparison['sample'] == metric['sample']
    assert comparison['phase'] == metric['phase']

def check_kinetic_changes(channel):
    density_history = channel['density_history']
    assert len(density_history) == channel['iterations_at_save']
    assert density_history[0]['availability'] == 'unavailable'
    for sample in density_history[1:]:
        assert sample['availability'] == 'measured'
        assert sample['value'] >= 0.0
        assert math.isclose(sample['effective_limit'], sample['absolute_limit'] + sample['relative_factor'] * sample['normalization_scale'], rel_tol=1e-15)
    history = channel['function_kinetic_history']
    first = history[0]['change_from_previous']
    assert first['availability'] == 'unavailable' and first['value'] is None
    for previous, current in zip(history, history[1:]):
        change = current['change_from_previous']
        now = sum(f['value'] for functions in current['blocks'].values() for f in functions)
        before = sum(f['value'] for functions in previous['blocks'].values() for f in functions)
        assert change['value'] == abs(now - before)
        assert current['trace']['value'] == now
        assert change['quantity'] == 'absolute_kinetic_trace_change'
        assert change['comparison'] == previous['sample']
        assert change['sample'] == current['sample']
        assert change['active_gate'] is False and change['target'] is None and change['satisfied'] is None
    final = channel['saved_state_validation']['measurements']['delta_kinetic_energy']
    assert final['valid_for_saved_state'] is True
    assert final['comparison'] == history[-2]['sample']
    assert abs(final['value'] - history[-1]['change_from_previous']['value']) < 1e-10
    # Step-input evidence must use its own preceding iterate, not the saved output.
    historical = channel['measurements']['delta_kinetic_energy']
    assert historical == history[-2]['change_from_previous']

fresh = work / 'Fresh'
stages = run('fresh', fresh)
# Exact reuse later overwrites Fresh/results.json with zero work counts. Keep
# the original fresh result for numerical comparisons and portable reporting.
(work / 'fresh-results.json').write_text(json.dumps({'stages': stages}, indent=2) + '\n')
gs, response = stages
assert 'convergence_report' in gs, 'ground convergence measurement report missing'
g = gs['convergence_report']
assert set(g['measurements']) == {'density_change', 'charge_density_change', 'spin_density_change', 'bsh_residual', 'energy_change', 'absolute_kinetic_trace_change'}
check_function_bsh(g['measurements']['bsh_residual'], 'max')
validation = g['saved_state_validation']
check_function_bsh(validation['measurements']['bsh_residual'], 'max')
assert gs['convergence_report']['published_in']['generation'] == gs['state_id']
assert validation['affects_stopping'] is False
for name, m in validation['measurements'].items():
    assert m['availability'] == 'measured', (name, m)
    assert m['valid_for_saved_state'] is True
    assert m['value'] >= 0
assert validation['comparison']['iteration'] < g['measurements']['energy_change']['sample']['iteration']
assert g['decided_by'] == 'engine_check'
assert 'state_id' not in g or g['state_id'] == gs['state_id']
for name, m in g['measurements'].items():
    assert m['availability'] == 'measured', (name, m)
    assert m['value'] >= 0 and m['phase']
    assert m['sample']['iteration'] >= 0
    assert m['valid_for_saved_state'] is False
assert g['measurements']['energy_change']['active_gate'] is False
assert g['measurements']['energy_change']['target'] == 1e-6
assert g['measurements']['absolute_kinetic_trace_change']['active_gate'] is False
r = response['convergence_report']
assert r['reference']['state_id'] == gs['state_id']
assert r['reference']['energy_change'] == g['measurements']['energy_change']
# The evidence names the reference generation it came from; the manifest's
# report is that generation's own (D-105).
assert {k: v for k, v in r['reference']['saved_state_validation'].items() if k != 'state_id'} == \
    {k: v for k, v in g['saved_state_validation'].items() if k != 'state_id'}
assert len(r['channels']) == 4  # static + dynamic, two supplied rungs
for c in r['channels']:
    assert 'state_id' not in c and c['archive']  # the generation names itself (D-105)
    m = c['measurements']
    assert m['density_change']['availability'] == 'measured'
    assert m['bsh_residual']['availability'] == 'measured'
    assert m['density_change']['effective_limit'] > 0
    assert m['bsh_residual']['valid_for_saved_state'] is False
    assert m['energy_change']['availability'] == 'not_applicable'
    assert m['delta_kinetic_energy']['availability'] == 'measured'
    check_function_bsh(m['bsh_residual'], 'quadrature')
    check_bsh_comparison(m['bsh_residual'])
    check_kinetic_changes(c)
    history = c['function_bsh_history']
    assert len(history) == c['iterations_at_save']
    assert [sample['sample']['iteration'] for sample in history] == list(range(c['iterations_at_save']))
    for sample in history:
        check_function_bsh(sample, 'quadrature')
        check_bsh_comparison(sample)
        assert sample['normalization_phase'] == 'solver_step_output_after_acceleration'
        assert sample['valid_for_saved_state'] is False
    check_function_bsh(c['saved_state_validation']['measurements']['bsh_residual'], 'quadrature')
    check_bsh_comparison(c['saved_state_validation']['measurements']['bsh_residual'])
    kinetics = c['function_kinetic_history']
    assert len(kinetics) == c['iterations_at_save'] + 1
    assert [k['sample']['iteration'] for k in kinetics] == list(range(c['iterations_at_save'] + 1))
    blocks = {'x_alpha'} if c['frequency'] == 0 else {'x_alpha', 'y_alpha'}
    for sample in history:
        assert set(sample['per_function']['blocks']) == blocks
    assert set(c['saved_state_validation']['measurements']['bsh_residual']['per_function']['blocks']) == blocks
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
    assert 'state_id' not in validation
    assert validation['affects_stopping'] is False
    for name in ('density_change', 'bsh_residual'):
        v = validation['measurements'][name]
        assert v['availability'] == 'measured' and v['valid_for_saved_state'] is True
        assert v['value'] >= 0 and v['effective_limit'] > 0
    assert validation['measurements']['delta_kinetic_energy']['availability'] == 'measured'
# First differences at every SCF rung must be unavailable, never invented zeros.
first = [row for row in gs_history(fresh) if row['kinetic_trace']['sample']['iteration'] == 0]
assert len(first) == 2
for row in first:
    evidence = row['measurements']
    assert evidence['bsh_residual']['availability'] == 'unavailable'
    assert 'per_function' not in evidence['bsh_residual']
    for name in ('energy_change', 'density_change', 'charge_density_change', 'spin_density_change', 'absolute_kinetic_trace_change'):
        assert evidence[name]['availability'] == 'unavailable' and evidence[name]['value'] is None
# Reported changes must agree with independently differenced recorded samples.
by_rung = {}
for row in gs_history(fresh):
    by_rung.setdefault(row['kinetic_trace']['sample']['threshold'], []).append(row)
for samples in by_rung.values():
    assert [r['kinetic_trace']['sample']['iteration'] for r in samples] == list(range(len(samples)))
    for before, after in zip(samples, samples[1:]):
        values = after['measurements']
        check_function_bsh(values['bsh_residual'], 'max')
        assert values['bsh_residual']['sample']['iteration'] == before['kinetic_trace']['sample']['iteration']
        delta = abs(after['energy']['value'] - before['energy']['value'])
        assert abs(values['energy_change']['value'] - delta) < 1e-14
        delta_t = abs(after['kinetic_trace']['value'] - before['kinetic_trace']['value'])
        assert abs(values['absolute_kinetic_trace_change']['value'] - delta_t) < 1e-14
# Final energy and kinetic differences must use the recorded prior iteration,
# never the almost-identical state just before final diagonalization.
v = g['saved_state_validation']
comparison = v['comparison']
previous = next(row for row in gs_history(fresh)
                if row['kinetic_trace']['sample']['threshold'] == comparison['threshold']
                and row['kinetic_trace']['sample']['iteration'] == comparison['iteration'])
assert comparison['kinetic_trace']['blocks'] == previous['kinetic_trace']['blocks']
assert comparison['energy'] == previous['energy']['value']
assert comparison['kinetic_energy'] == previous['occupation_weighted_physical_kinetic_energy']['value']
for metric, key in [('energy_change', 'energy')]:
    assert abs(v['measurements'][metric]['value'] - abs(v[key] - comparison[key])) < 1e-14
assert abs(v['measurements']['absolute_kinetic_trace_change']['value'] -
           abs(v['kinetic_trace']['value'] - comparison['kinetic_trace']['value'])) < 1e-14
for sample in g['iteration_history']:
    trace = sample['kinetic_trace']
    assert trace['value'] == sum(f['value'] for block in trace['blocks'].values() for f in block)
    assert trace['weighting'] == 'none'
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
for channel in resumed['convergence_report']['channels']:
    check_kinetic_changes(channel)
    history = channel['function_bsh_history']
    assert [sample['sample']['iteration'] for sample in history] == list(range(channel['iterations_at_save']))
    for sample in history:
        check_function_bsh(sample, 'quadrature')
        check_bsh_comparison(sample)
        assert sample['normalization_phase'] == 'solver_step_output_after_acceleration'
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
    a = one['saved_state_validation']['measurements']['bsh_residual']['per_function']['blocks']
    b = two['measurements']['bsh_residual']['per_function']['blocks']
    for block in a:
        for x, y in zip(a[block], b[block]):
            assert x['function_index'] == y['function_index']
            assert abs(x['value'] - y['value']) < 1e-10
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
for x, y in zip(one['bsh_residual']['per_function']['blocks']['alpha'],
                three['bsh_residual']['per_function']['blocks']['alpha']):
    assert x['function_index'] == y['function_index']
    assert abs(x['value'] - y['value']) < 1e-9
step2 = next(row for row in gs_history(work / 'GroundSteps3')
             if row['kinetic_trace']['sample']['iteration'] == 1)
for metric in ('density_change', 'charge_density_change', 'spin_density_change', 'energy_change', 'absolute_kinetic_trace_change'):
    expected = step2['measurements'][metric]['value']
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
# Trace aggregates and temporal changes are scientific evidence, not trusted flags.
for location in ('raw_trace', 'trace_change', 'history_trace', 'input_phase', 'history_phase', 'saved_phase'):
    root = work / ('InvalidTrace-' + location)
    shutil.copytree(fresh, root)
    store = root / 'gs/checkpoints'
    latest = (store / 'latest').read_text()
    mfile = store / json.loads(latest)['generation'] / 'manifest.json'
    m = json.loads(mfile.read_text())
    report = m['metrics']['convergence_report']
    if location == 'input_phase':
        m['metrics']['kinetic_trace']['phase'] = 'incompatible_phase'
    elif location == 'history_phase':
        report['iteration_history'][0]['kinetic_trace']['phase'] = 'incompatible_phase'
    elif location == 'saved_phase':
        report['saved_state_validation']['kinetic_trace']['phase'] = 'incompatible_phase'
    elif location == 'raw_trace':
        m['metrics']['kinetic_trace']['value'] = 99.0
    elif location == 'history_trace':
        report['iteration_history'][0]['kinetic_trace']['value'] = 99.0
    else:
        report['saved_state_validation']['measurements']['absolute_kinetic_trace_change']['value'] = 99.0
    mfile.write_text(json.dumps(m))
    run('invalid-trace-' + location, root, resume=True, success=False)
    assert 'invalid ground-state measurement' in (work / ('invalid-trace-' + location + '.log')).read_text()
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
# Derived diagnostics must be checked against their actual source samples.
for corrupt in ('maximum', 'kinetic_change', 'kinetic_sample', 'historical_kinetic', 'bsh_limit'):
    root = work / ('InvalidComparison-' + corrupt)
    shutil.copytree(fresh, root)
    store = root / 'optical/checkpoints'
    latest = (store / 'latest').read_text()
    mfile = store / json.loads(latest)['generation'] / 'manifest.json'
    manifest = json.loads(mfile.read_text())
    channels = manifest['numerical']['engine_metadata']['fd_states']
    entry = next(iter(next(iter(next(iter(channels.values())).values())).values()))
    if corrupt == 'maximum':
        entry['function_bsh_history'][0]['maximum_comparison']['value'] = 0.0
    elif corrupt == 'kinetic_change':
        entry['function_kinetic_history'][1]['change_from_previous']['current_trace'] = 99.0
    elif corrupt == 'historical_kinetic':
        entry['convergence_measurements']['kinetic_energy']['blocks']['x_alpha'][0]['value'] = 77.0
    elif corrupt == 'kinetic_sample':
        entry['function_kinetic_history'][1]['sample']['k'] += 2
    else:
        entry['function_bsh_history'][0]['effective_limit'] = 99.0
    mfile.write_text(json.dumps(manifest))
    run('invalid-comparison-' + corrupt, root, resume=True, success=False)
    # The response dconv is a written target (DESIGN B8), so the written-target
    # validator checks every named response entry; it runs before the
    # metric-comparison check and is the one that rejects a changed sample k.
    expected_error = {'historical_kinetic': 'per-function kinetic measurement',
                      'kinetic_sample': 'explicit accuracy evidence'}.get(corrupt, 'metric comparison')
    assert 'invalid response state: ' + expected_error in (work / ('invalid-comparison-' + corrupt + '.log')).read_text()
    assert (store / 'latest').read_text() == latest

# Reject corrupt component values in both historical and saved-state evidence.
for stage in ('gs', 'optical'):
    for location in ('historical', 'saved'):
        root = work / ('InvalidFunctionBSH-' + stage + '-' + location)
        shutil.copytree(fresh, root)
        store = root / stage / 'checkpoints'
        latest = (store / 'latest').read_text()
        mfile = store / json.loads(latest)['generation'] / 'manifest.json'
        manifest = json.loads(mfile.read_text())
        if stage == 'gs':
            entry = manifest['metrics']['convergence_report']
            measurements = entry['measurements']
        else:
            channels = manifest['numerical']['engine_metadata']['fd_states']
            entry = next(iter(next(iter(next(iter(channels.values())).values())).values()))
            measurements = entry['convergence_measurements']
        if location == 'saved':
            measurements = entry['saved_state_validation']['measurements']
        blocks = measurements['bsh_residual']['per_function']['blocks']
        next(iter(blocks.values()))[0]['value'] = None
        mfile.write_text(json.dumps(manifest))
        label = 'invalid-function-bsh-' + stage + '-' + location
        run(label, root, resume=True, success=False)
        assert 'per-function BSH' in (work / (label + '.log')).read_text()
        assert (store / 'latest').read_text() == latest

# History contains its own aggregate as well as components; validate both.
for corrupt in ('aggregate', 'component'):
    root = work / ('InvalidBSHHistory-' + corrupt)
    shutil.copytree(fresh, root)
    store = root / 'optical/checkpoints'
    latest = (store / 'latest').read_text()
    mfile = store / json.loads(latest)['generation'] / 'manifest.json'
    manifest = json.loads(mfile.read_text())
    channels = manifest['numerical']['engine_metadata']['fd_states']
    entry = next(iter(next(iter(next(iter(channels.values())).values())).values()))
    sample = entry['function_bsh_history'][0]
    if corrupt == 'aggregate':
        sample['value'] = None
    else:
        sample['per_function']['blocks']['x_alpha'][0]['value'] = None
    mfile.write_text(json.dumps(manifest))
    label = 'invalid-bsh-history-' + corrupt
    run(label, root, resume=True, success=False)
    assert 'per-function BSH' in (work / (label + '.log')).read_text()
    assert (store / 'latest').read_text() == latest

# Aggregates from before this extension remain usable without invented components.
legacy = work / 'LegacyAggregates'
shutil.copytree(fresh, legacy)
for stage in ('gs', 'optical'):
    store = legacy / stage / 'checkpoints'
    mfile = store / json.loads((store / 'latest').read_text())['generation'] / 'manifest.json'
    manifest = json.loads(mfile.read_text())
    if stage == 'gs':
        # Written before S4a too: no state record (numerical.states).
        manifest['numerical'].pop('states')
        report = manifest['metrics']['convergence_report']
        report['measurements']['bsh_residual'].pop('per_function')
        report['saved_state_validation']['measurements']['bsh_residual'].pop('per_function')
    else:
        for channel in manifest['numerical']['engine_metadata']['fd_states'].values():
            for rung in channel.values():
                for entry in rung.values():
                    strip_written_targets(entry)
                    entry['convergence_measurements']['bsh_residual'].pop('per_function')
                    entry['convergence_measurements']['bsh_residual'].pop('maximum_comparison', None)
                    entry['saved_state_validation']['measurements']['bsh_residual'].pop('per_function')
                    entry['saved_state_validation']['measurements']['bsh_residual'].pop('maximum_comparison', None)
                    entry.pop('function_bsh_history')
        # A response generation from then has no published report (D-105).
        manifest['metrics'].pop('convergence_report', None)
    mfile.write_text(json.dumps(manifest))
legacy_results = run('legacy-aggregates', legacy, resume=True)
assert all(s['iterations_this_run'] == 0 for s in legacy_results)
legacy_ground = legacy_results[0]['convergence_report']
assert legacy_ground['measurements']['bsh_residual']['availability'] == 'measured'
assert 'per_function' not in legacy_ground['measurements']['bsh_residual']
assert 'per_function' not in legacy_ground['saved_state_validation']['measurements']['bsh_residual']
for channel in legacy_results[1]['convergence_report']['channels']:
    assert channel['measurements']['bsh_residual']['availability'] == 'measured'
    assert 'per_function' not in channel['measurements']['bsh_residual']
    assert 'per_function' not in channel['saved_state_validation']['measurements']['bsh_residual']
    assert 'function_bsh_history' not in channel

# Reporting-only generations from before metric comparisons retain their raw
# histories, without acquiring fabricated comparisons during exact reuse.
legacy = work / 'LegacyComparisons'
shutil.copytree(fresh, legacy)
store = legacy / 'optical/checkpoints'
mfile = store / json.loads((store / 'latest').read_text())['generation'] / 'manifest.json'
manifest = json.loads(mfile.read_text())
for channel in manifest['numerical']['engine_metadata']['fd_states'].values():
    for rung in channel.values():
        for entry in rung.values():
            strip_written_targets(entry)
            for measurements in (entry['convergence_measurements'], entry['saved_state_validation']['measurements']):
                measurements['bsh_residual'].pop('maximum_comparison')
                measurements['delta_kinetic_energy'] = dict(value=None, availability='not_applicable', active_gate=False,
                    reason='kinetic_values_reported_per_function_without_scalar_reduction')
            for sample in entry['function_kinetic_history']:
                sample.pop('change_from_previous')
            for sample in entry['function_bsh_history']:
                for key in ('maximum_comparison', 'absolute_limit', 'relative_factor', 'normalization_scale',
                            'normalization_phase', 'effective_limit', 'satisfied'):
                    sample.pop(key)
# A response generation from then has no published report (D-105).
manifest['metrics'].pop('convergence_report', None)
mfile.write_text(json.dumps(manifest))
legacy_results = run('legacy-comparisons', legacy, resume=True)
assert all(stage['iterations_this_run'] == 0 for stage in legacy_results)
for channel in legacy_results[1]['convergence_report']['channels']:
    assert 'maximum_comparison' not in channel['measurements']['bsh_residual']
    assert all('change_from_previous' not in sample for sample in channel['function_kinetic_history'])
    assert all('effective_limit' not in sample for sample in channel['function_bsh_history'])

# Old generations stay reusable but cannot acquire invented measurements.
old = work / 'Old'
shutil.copytree(fresh, old)
for stage in ('gs', 'optical'):
    store = old / stage / 'checkpoints'
    mfile = store / json.loads((store / 'latest').read_text())['generation'] / 'manifest.json'
    m = json.loads(mfile.read_text())
    if stage == 'gs':
        # Written before S4a too: no state record (numerical.states).
        m['numerical'].pop('states')
        m['metrics'].pop('convergence_report')
    else:
        for channel in m['numerical']['engine_metadata']['fd_states'].values():
            for rung in channel.values():
                for e in rung.values():
                    strip_written_targets(e)
                    e.pop('convergence_measurements')
                    e.pop('saved_state_validation', None)
                    e.pop('function_kinetic_history', None)
                    e.pop('function_bsh_history', None)
        # Nor did it publish a report (D-105).
        m['metrics'].pop('convergence_report', None)
    mfile.write_text(json.dumps(m))
# Optional legacy reports remain optional only when no explicit policy needs them.
for control in ('accuracy', 'progress', 'progress-only'):
    required = work / ('RequiredEvidence-' + control)
    shutil.copytree(old, required)
    store = required / 'gs/checkpoints'
    pointer = (store / 'latest').read_text()
    mfile = store / json.loads(pointer)['generation'] / 'manifest.json'
    manifest = json.loads(mfile.read_text())
    controls = manifest['numerical']['parameters']['parameters']
    if control != 'progress-only':
        controls['accuracy.kinetic_trace_change'] = 0.0
    if control.startswith('progress'):
        controls['progress.window'] = 1
        controls['progress.minimum_relative_improvement'] = 0.1
    mfile.write_text(json.dumps(manifest))
    run('required-evidence-' + control, required, resume=True, success=False)
    assert (store / 'latest').read_text() == pointer
    assert 'missing configured accuracy/progress report' in (work / ('required-evidence-' + control + '.log')).read_text()

# A response generation from before these reports: still reused, and no
# measurement is invented for it (its reference ground state is current).
old_response = work / 'OldResponse'
shutil.copytree(fresh, old_response)
shutil.rmtree(old_response / 'optical')
shutil.copytree(old / 'optical', old_response / 'optical')
old_results = run('old-response', old_response, resume=True)
assert all(s['iterations_this_run'] == 0 for s in old_results)
assert old_results[1]['convergence_report']['channels'][0]['measurements']['bsh_residual']['value'] is None
assert all('saved_state_validation' not in c and 'function_kinetic_history' not in c and 'function_bsh_history' not in c
           for c in old_results[1]['convergence_report']['channels'])
# A ground state from before the report has no measured values to re-judge
# against today's request, so --continue resumes it rather than trusting its
# old verdict (DESIGN B6: never a saved verdict; the moldft counterpart of
# D-005.4). Its new generations carry the measurements.
old_ground = work / 'OldGround'
shutil.copytree(old, old_ground)
shutil.rmtree(old_ground / 'optical')
ground = run('old-ground', old_ground, deck.split('molresponse optical')[0], resume=True)[0]
assert ground['work'] == 'resumed' and ground['iterations_this_run'] > 0, ground['work']
assert 'cannot be re-judged' in ground['continuation']['why'], ground['continuation']
assert ground['convergence_report']['measurements']['energy_change']['availability'] == 'measured'
print('named convergence reporting, saved-state validation and continuation passed', flush=True)
