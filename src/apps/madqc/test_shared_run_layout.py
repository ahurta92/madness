#!/usr/bin/env python3
"""Run-root records, byte-preserving offline relocation and native restart."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
work = Path(tempfile.mkdtemp(prefix='shared-run-layout-'))
print('artifacts:', work, flush=True)
(work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 l 30
 plotdens true
 plotcube true
 npt_plot 5
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 0
end
response
 dipole true
 dipole.directions z
 dipole.frequencies [0.0]
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 0
end
molecule
 units atomic
 He 0 0 0
end
''')
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
root = work/'nested/run'
def run(name, extra=(), success=True):
    result = subprocess.run(launch+[str(Path(a.binary).resolve()), 'input.in', '--wf=response',
                            '--prefix=nested/run', *extra], cwd=work, env=env,
                            capture_output=True, text=True, timeout=180)
    output = result.stdout + result.stderr
    (work/(name+'.log')).write_text(output)
    assert (result.returncode == 0) == success, output[-3000:]
    return output

def hashes(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob('*') if p.is_file()}

def native():
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (root/'checkpoints').rglob('*')
            if p.is_file() and (p.name == 'latest' or p.suffix == '.00000' or p.name == 'manifest.json')}

run('fresh')
assert (root/'results.json').is_file(), 'shared results must live at the run root'
assert (root/'timings.json').is_file()
assert (root/'checkpoints/run.json').is_file()
assert not list(root.glob('task_*'))
assert (root/'analysis/ground_state/moldft.calc_info.json').is_file()
assert (root/'analysis/response/response_metadata.json').is_file()
assert list((root/'analysis/response').glob('*.jsonl'))
assert not list((root/'analysis').rglob('checkpoints'))
viz = json.loads((work/'nested/run.viz_manifest.json').read_text())
plots = [p for p in viz['artifacts'] if p['format'] == 'cube']
assert plots and all(p['calculation_id'] == 'ground_state' for p in plots)
assert all((work/p['path']).is_file() for p in plots)
original = native()
source_records = {p.relative_to(root): p.read_bytes()
                  for p in (root/'checkpoints').rglob('result_sources.json')}
run('continue', ['--continue'])
assert native() == original, 'reuse should preserve the exact native sources'
# A marker also protects the interval after all moves but before cleanup.
marker = root/'.madqc-layout-migration.json'
marker.write_text('{}')
before = hashes(root)
assert 'migrate_madqc_run.py' in run('partial_move', ['--continue'], success=False)
assert hashes(root) == before
marker.unlink()

# Reconstruct the previous development layout without changing archived bytes.
old = root/'task_0/moldft'
(root/'task_0').mkdir()
(root/'analysis/ground_state').rename(old)
for name in ['checkpoints', 'results.json', 'timings.json']:
    (root/name).rename(old/name)
(root/'task_1').mkdir()
(root/'analysis/response').rename(root/'task_1/molresponse')
response = root/'task_1/molresponse'
(response/'analysis').mkdir()
for f in response.glob('*.jsonl'):
    f.rename(response/'analysis'/f.name)
root.joinpath('analysis').rmdir()
before = hashes(root)
for name, extra in [('old_auto', []), ('old_continue', ['--continue']),
                    ('old_none', ['--dft=restart=none']),
                    ('old_inspect', ['--wf=scf', '--dft=restart=read_only'])]:
    output = run(name, extra, success=False)
    assert 'migrate_madqc_run.py' in output
    assert hashes(root) == before, 'old-layout rejection must not alter saved work'

migration = Path(__file__).resolve().parents[3]/'bin/migrate_madqc_run.py'
result = subprocess.run([sys.executable, str(migration), str(root)], capture_output=True, text=True)
assert result.returncode == 0, result.stdout + result.stderr
assert native() == original
assert not list(root.glob('task_*'))
assert not (root/'.madqc-layout-migration.json').exists()
run('relocated_continue', ['--continue'])
assert native() == original, 'relocated checkpoints must reuse without recomputation'
assert all((root/p).read_bytes() == value for p, value in source_records.items()), 'old results lost their original source records'
# Inspection uses the root and leaves every saved record untouched.
before = hashes(root)
run('inspect', ['--wf=scf', '--dft=restart=read_only'])
assert hashes(root) == before
print('shared run layout checks passed', flush=True)
