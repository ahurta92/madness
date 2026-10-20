#!/usr/bin/env python3
"""Catch missing publication, false convergence, mixed generations and unsafe reuse."""
import argparse
import hashlib
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
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='named-checkpoint-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
deck = (Path(a.source) / 'named_h2.in').read_text().split('molresponse optical')[0]

def run(label, root, text=deck, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=240)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()

def published(root):
    store = root / 'gs/checkpoints'
    assert (store / 'latest').exists(), 'moldft did not publish a checkpoint'
    generation = json.loads((store / 'latest').read_text())['generation']
    manifest = json.loads((store / generation / 'manifest.json').read_text())
    return store / generation, manifest

def result(root):
    return json.loads((root / 'gs/results.json').read_text())

def hashes(directory):
    return {str(f.relative_to(directory)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in directory.rglob('*') if f.is_file()}

fresh = work / 'Fresh'
run('fresh', fresh)
generation, manifest = published(fresh)
assert len(manifest['archives']) == 1 and manifest['archives'][0]['role'] == 'orbitals'
assert result(fresh)['checkpoint_generation'] == manifest['generation']
assert manifest['stop_reason'] == 'converged'
assert all(v['generation'] == manifest['generation'] for v in manifest['properties'])
original = hashes(generation)
fresh_iterations = result(fresh)['iterations_this_run']
energy = result(fresh)['properties']['energy']['value']
assert abs(energy - (-1.1336281435304234)) < 2e-7
assert max(abs(x) for x in result(fresh)['properties']['dipole']['value']) < 1e-8
run('reuse', fresh, resume=True)
assert result(fresh)['work'] == 'loaded'
assert result(fresh)['iterations_this_run'] == 0
assert published(fresh)[0] == generation and hashes(generation) == original

# Changing a target forces iteration; old convergence cannot certify an update
# that never reached the solver's convergence test (which starts at iteration 1).
stale = work / 'StaleEvidence'
shutil.copytree(fresh, stale)
run('stale-evidence', stale, deck.replace('econv 1e-6', 'econv 1e-7').replace('maxiter 30', 'maxiter 1'),
    resume=True, success=False)
assert result(stale)['status'] == 'unconverged'
assert published(stale)[1]['stop_reason'] != 'converged'

# Unpublished files must never be selected, even when newer than latest.
(fresh / 'gs/checkpoints/incomplete').mkdir()
(fresh / 'gs/checkpoints/incomplete/state.restartdata.00000').write_bytes(b'incomplete')
run('ignore-unpublished', fresh, resume=True)
assert published(fresh)[0] == generation

# A failed publisher may leave an unpublished generation, but cannot replace latest.
lock = fresh / 'gs/checkpoints/.publish-lock'
lock.mkdir()
log = run('publication-failure', fresh, deck.replace('dconv 1e-5', 'dconv 1e-6'),
          resume=True, success=False)
assert '.publish-lock' in log
assert published(fresh)[0] == generation and hashes(generation) == original
lock.rmdir()
run('reuse-after-failure', fresh, resume=True)
assert published(fresh)[0] == generation

partial = work / 'Partial'
run('limited', partial, (Path(a.source) / 'named_h2.in').read_text().replace('maxiter 30', 'maxiter 2'), success=False)
old_dir, old = published(partial)
assert old['stop_reason'] != 'converged'
assert result(partial)['status'] == 'unconverged'
assert not list((partial / 'optical').iterdir()), 'response ran against an unfinished reference'
old_hash = hashes(old_dir)
run('resume', partial, resume=True)
assert result(partial)['work'] == 'resumed'
assert result(partial)['acceleration_history'] == 'reset'
assert abs(result(partial)['properties']['energy']['value'] - energy) < 2e-7
assert hashes(old_dir) == old_hash

coarse = work / 'Coarse'
run('coarse', coarse, deck.replace('[1e-4, 1e-6]', '[1e-4]').replace('dconv 1e-5', 'dconv 1e-3'))
coarse_dir, _ = published(coarse)
coarse_hash = hashes(coarse_dir)
run('tighten', coarse, resume=True)
assert result(coarse)['work'] == 'resumed'
assert abs(result(coarse)['properties']['energy']['value'] - energy) < 2e-7
assert hashes(coarse_dir) == coarse_hash

# Full-inventory validation must fail before any replacement generation.
for defect in ('missing', 'corrupt', 'inventory', 'geometry'):
    broken = work / defect
    shutil.copytree(fresh, broken)
    directory, m = published(broken)
    shard = directory / m['archives'][0]['shards'][-1]['path']
    if defect == 'missing':
        shard.unlink()
    elif defect == 'corrupt':
        with shard.open('r+b') as f:
            f.seek(-1, 2)
            b = f.read(1)
            f.seek(-1, 2)
            f.write(bytes([b[0] ^ 1]))
    elif defect == 'inventory':
        m['archives'] = []
        (directory / 'manifest.json').write_text(json.dumps(m))
    before = hashes(broken / 'gs/checkpoints')
    text = deck.replace('H 0 0 0.7', 'H 0 0 0.8') if defect == 'geometry' else deck
    log = run(defect, broken, text, resume=True, success=False)
    assert 'checkpoint' in log.lower(), log[-3000:]
    assert hashes(broken / 'gs/checkpoints') == before

print(json.dumps({'mpi_ranks': a.ranks, 'energy': energy,
                  'fresh_iterations': fresh_iterations,
                  'resumed_iterations': result(partial)['iterations_this_run'],
                  'tightened_iterations': result(coarse)['iterations_this_run']}, indent=2))
