#!/usr/bin/env python3
"""Single alpha state: saved-state evidence, reuse and tighter continuation.

Loose/zero gates distinguish lifecycle decisions. A separate final run checks
alpha against the unchanged helium reference; no reference regeneration.
"""
from response_test_checks import assert_response_setup_timed
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--frequency', type=float, default=0.05)
p.add_argument('--mpiexec', default='mpiexec')
p.add_argument('--wiring-only', action='store_true')
a = p.parse_args()
work = Path(tempfile.mkdtemp(prefix='shared-response-'))
print('artifacts:', work, flush=True)
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
binary = str(Path(a.binary).resolve())
root = work/'run'
task = f'fd:dipole_z@f{a.frequency:.5f}'
tasks = root/'checkpoints/response/tasks'/task
counter = 0


def run(settings, success=True, extra=()):
    global counter
    (work/'input.in').write_text('''dft
 xc hf
 localize canon
 k 6
 protocol [1e-4]
 maxiter 25
 convergence.bsh_residual 1e-4
 convergence.density_change 1e-4
 convergence.energy_change 1e-6
 print_level 0
end
response
 dipole.directions z
 dipole.frequencies [0.05]
 print_level 1
 kain true
 convergence.bsh_residual 1000
 convergence.density_change 1000
 maxiter 1
''' + settings + '''
end
molecule
 units atomic
 He 0 0 0
end
''')
    deck=work/'input.in'
    deck.write_text(deck.read_text().replace('dipole.frequencies [0.05]',f'dipole.frequencies [{a.frequency}]'))
    result = subprocess.run(launch+[binary,'input.in','--wf=response','--prefix=run']+list(extra),
        cwd=work,env=env,capture_output=True,text=True,timeout=240)
    counter += 1
    (work/f'{counter:02d}.log').write_text(result.stdout+result.stderr)
    assert (result.returncode == 0) == success, (result.returncode,result.stdout[-5000:],result.stderr)
    return result.stdout


def saved():
    generation = json.loads((tasks/'latest').read_text())['generation']
    return generation,json.loads((tasks/generation/'manifest.json').read_text())


def metric(record, name):
    # Shared manifest schema, independent of the legacy response analysis file.
    return next(m for m in record['task']['iteration']['metrics'] if m['metric']==name)


def response_summary():
    metadata = next((work/'run').rglob('response_metadata.json'))
    return json.loads(metadata.read_text())['fd_states']['dipole_z']['1e-04_k6'][f'f{a.frequency:.5f}']

def current_timing():
    records = json.loads((root/'timings.json').read_text())
    assert records['version'] == 1
    current = records['attempts'][-1]
    assert current['execution'] == {'mpi_ranks': a.ranks, 'threads_per_rank': 2}
    assert current['timing']['clock_scope'] == 'driver_rank'
    assert_response_setup_timed(current)
    def walk(node):
        assert node['name']
        if 'seconds' in node:
            assert node['seconds'] >= 0
        else:
            assert node['availability']['status'] == 'unavailable'
        yield node
        for child in node.get('children', []):
            yield from walk(child)
    return records, list(walk(current['timing']))

first_output = run('')
first_timing, first_steps = current_timing()
# One update with loose targets: change metrics need a predecessor.
assert [n['iteration'] for n in first_steps if n['name'] == 'iteration'] == [0, 1]
assert any(n['name'] == 'prepare' for n in first_steps)
assert any(n['name'] == 'checkpoint/report' for n in first_steps)
assert any(n['name'] == 'advance' for n in first_steps)
assert 'policy: dconv_user' not in first_output, 'shared response must report its actual criteria'
summary = first_output.split('MADQC RESULTS SUMMARY')[-1]
assert 'dconv =' not in summary and '1e+10' not in summary, 'legacy stamps leaked into summary'
assert 'Wall time' in summary and '2 threads' in summary
assert 'Response steps' in summary and 'solve ' in summary
assert 'Timing response' in first_output and 'Evaluation' in first_output and 'Checkpoint and reporting' in first_output
assert first_output.count('Solver update') == 1, 'one update reports one update timing'
assert 'Induced-density change' in first_output and 'Response kinetic change' in first_output
assert 'Checkpoint saved.' in first_output
report = json.loads((work/'run.calc_info.json').read_text())
for application in report['tasks']:
    runtime = application['metadata']
    assert runtime['elapsed_time'] > 0 and runtime['clock_scope'] == 'driver_rank'
    assert runtime['mpi_size'] == a.ranks and runtime['nthreads'] == 2

first, record = saved()
assert record['task']['iteration']['updates_this_attempt'] == 1
assert metric(record,'density_change')['status'] == 'measured'
assert metric(record,'energy_change')['status'] == 'not_applicable'
initial_residual = metric(record,'bsh_residual')['value']
assert initial_residual > 0
original = (tasks/first/'manifest.json').read_bytes()
archives = {str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (tasks/first).glob('*.00000')}
assert archives
reuse_output = run('', extra=['--response=convergence.bsh_residual=500'])
reuse_timing, reuse_steps = current_timing()
assert reuse_timing['attempts'][:len(first_timing['attempts'])] == first_timing['attempts']
new_attempts = reuse_timing['attempts'][len(first_timing['attempts']):]
assert len(new_attempts) == 2
assert {a['timing']['children'][0]['calculation_id'] for a in new_attempts} == {'ground_state', 'response'}
assert not any(n['name'] == 'iteration' for n in reuse_steps)
assert any(n['name'] == 'load' for n in reuse_steps)
assert 'Timing response' not in reuse_output, 'reuse must not replay historical solver timings'
assert 'Wall time' in reuse_output
assert 'Checkpoint response' in reuse_output and 'Load and checks' in reuse_output
assert 'Saved metrics against current tolerances' in reuse_output
assert 'no solver iterations required' in reuse_output
assert 'Checkpoint saved.' not in reuse_output
assert saved()[0] == first, 'satisfying tighter tolerance must reuse saved evidence'
assert response_summary()['shared_evaluations'] == 0
assert response_summary()['shared_reused']
assert all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in archives.items())
run('convergence.bsh_residual 0\nmaxiter 0', success=False)
second, resumed = saved()
assert second != first
assert resumed['task']['iteration']['updates_this_attempt'] == 0
assert abs(metric(resumed,'bsh_residual')['value']-initial_residual) < 1e-10
assert metric(resumed,'density_change')['status'] == 'unavailable'
update_output = run('convergence.bsh_residual 0\nmaxiter 1', success=False)
assert update_output.count('Timing response') == 2
assert update_output.count('Solver update') == 1
third, updated = saved()
assert third != second and updated['task']['iteration']['updates_this_attempt'] == 1
assert metric(updated,'density_change')['status'] == 'measured'
updated_residual = metric(updated,'bsh_residual')['value']
assert updated_residual != initial_residual
run('convergence.bsh_residual 0\nmaxiter 0', success=False)
fourth, checked = saved()
assert fourth != third
assert abs(metric(checked,'bsh_residual')['value']-updated_residual) < 1e-10, 'checkpoint metrics describe different functions'
assert metric(checked,'density_change')['status'] == 'unavailable'
assert abs(metric(checked,'kinetic_change')['evidence']['current']-
           metric(updated,'kinetic_change')['evidence']['current']) < 1e-10
run('convergence.bsh_residual 0\nmaxiter 2\ncheckpoint.every_iterations 1', success=False)
assert saved()[1]['task']['iteration']['updates_this_attempt'] == 2
assert len(list(tasks.glob('*/*.00000'))) == 1, 'intermediate payloads must be reclaimed; audit manifests remain'
if a.wiring_only:
    print('response checkpoint/result/timing wiring passed', flush=True)
    raise SystemExit(0)
run('convergence.bsh_residual 1e-5\nconvergence.density_change 1e-5\nconvergence.kinetic_change 1e-5\nmaxiter 30')
final, final_record = saved()
assert final_record['task']['outcome']['complete']
assert metric(final_record,'density_change')['status'] == 'measured'
assert metric(final_record,'kinetic_change')['status'] == 'measured'
results = json.loads((root/'results.json').read_text())['calculations']
assert {x['calculation_id'] for x in results} == {'ground_state','response'}
response = next(x for x in results if x['calculation_id']=='response')
alpha = response['properties']['polarizability']['values'][0]['components']['zz']
reference = json.loads((Path(a.source)/'src/examples/qc/response_he_alpha/reference/response_he_alpha.calc_info.json').read_text())
expected = reference['tasks'][1]['properties']['response_properties']['alpha']['1e-04_k6'][0 if a.frequency==0 else 1]['alpha'][0][0]
assert abs(alpha-expected) < 1e-3, (alpha,expected)
print('alpha:', alpha, 'reference:', expected, 'absolute difference:', abs(alpha-expected), flush=True)
assert final_record['task']['inputs'][0]['calculation_id'] == 'ground_state'
assert (tasks/first/'.').exists() is False or (tasks/first/'manifest.json').read_bytes()==original
for history in (work/'run/analysis/response').glob('*.jsonl'):
    history.unlink()
if a.wiring_only:
    print('response checkpoint/result/timing wiring passed', flush=True)
    raise SystemExit(0)
run('convergence.bsh_residual 1e-5\nconvergence.density_change 1e-5\nconvergence.kinetic_change 1e-5\nmaxiter 30')
assert saved()[0] == final
assert response_summary()['shared_evaluations'] == 0
# Retired keys fail with a migration message and leave the saved state alone.
for retired, message in [('convergence.criteria [bsh_residual]', 'convergence.criteria has been retired'),
                         ('convergence.bsh_residual.relative 1', 'convergence.bsh_residual.relative has been retired')]:
    run(retired, success=False)
    assert message in (work/f'{counter:02d}.log').read_text(), retired
    assert saved()[0] == final
for invalid in ['convergence.bsh_residual -1\nconvergence.bsh_residual 1000',
                'excited.maxiter 1', 'dipole.directions xq']:
    run(invalid, success=False)
    assert saved()[0] == final
run('', success=False, extra=['--response=convergence.bsh_residual=-1',
                             '--response=convergence.bsh_residual=1000'])
run('', success=False, extra=['--response=unknown_criterion=1', '--response=maxiter=0'])
assert saved()[0] == final
payload = next((tasks/final).glob('*.00000'))
with payload.open('r+b') as stream:
    stream.seek(-1,2); value=stream.read(1); stream.seek(-1,2); stream.write(bytes([value[0]^1]))
run('convergence.bsh_residual 1000', success=False)
assert saved()[0] == final
print('shared driven alpha and tighter restart passed')
