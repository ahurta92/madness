#!/usr/bin/env python3
"""A moldft stage seeded from a Dalton directory (D-112): its orbitals start
from the Dalton molden.inp when no saved state is used, as the unnamed path's
ground-state seed does. The twin: water at the geometry of the qctest
response_h2o_tpa_seeded, seeded from that case's dalton/ and cold, the same
answer in fewer SCF iterations. Then: every generation holds the seed files,
listed in the manifest's seed block; a --continue loads the saved state and
seeds nothing; an altered seed copy is refused on read; another geometry is
refused before the SCF."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)
import named_workdir  # noqa: E402  (D-108: removed when the test passes)

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = named_workdir.work_directory('named-moldft-seed-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
dalton = (Path(a.source).parent.parent / 'examples/qc/response_h2o_tpa_seeded/dalton').resolve()
assert (dalton / 'molden.inp').exists(), dalton

molecule = """molecule water
  eprec 1e-6
  units atomic
  no_orient true
  O   0.0000000000   0.0000005943   0.1925923500
  H   0.0000000000   1.4214240411  -0.8721230117
  H   0.0000000000  -1.4214246355  -0.8721234678
end
"""
cold = molecule + """
moldft gs
  molecule water
  xc hf
  L 200
  localize new
  protocol [1e-4]
  dconv 1e-4
  maxiter 30
  requested_properties [energy]
end
"""
seeded = cold.replace('  maxiter 30\n', '  maxiter 30\n  seed.dalton_dir ' + str(dalton) + '\n')
assert seeded != cold


def run(label, root, text, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=900)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


def result(root):
    return named_records.stage_record(root / 'gs')


def manifest(root):
    store = root / 'gs/checkpoints'
    generation = json.loads((store / 'latest').read_text())['generation']
    return generation, json.loads((store / generation / 'manifest.json').read_text())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# The twin: the same answer, fewer SCF iterations.
run('cold', work / 'Cold', cold)
log = run('seeded', work / 'Seeded', seeded)
c, s = result(work / 'Cold'), result(work / 'Seeded')
assert c['status'] == s['status'] == 'completed', (c['status'], s['status'])
cold_iterations, seeded_iterations = c['iterations_this_run'], s['iterations_this_run']
print('SCF iterations: cold', cold_iterations, 'seeded', seeded_iterations, flush=True)
assert seeded_iterations < cold_iterations, (seeded_iterations, cold_iterations)
energy_c, energy_s = c['properties']['energy']['value'], s['properties']['energy']['value']
print('energies: cold', energy_c, 'seeded', energy_s, 'difference', energy_s - energy_c, flush=True)
assert abs(energy_s - energy_c) < 5e-5, (energy_c, energy_s)
assert '[DALTON-SEED] GS: seeding' in log, log[-3000:]
# No .out in the qctest's dalton/: the method could not be checked, and the
# stage says so (D-112 condition 5).
assert 'stage gs: Dalton seed: no Dalton .out in the directory' in log, log[-3000:]

# The record (no publish-sample case holds a seeded moldft stage, so the shape
# is pinned here): every generation holds seed/molden.inp, listed by size and
# SHA-256, with what the seed recorded.
generation, m = manifest(work / 'Seeded')
seed = m['seed']
assert seed['dalton_dir'] == str(dalton), seed['dalton_dir']
assert [f['path'] for f in seed['files']] == ['seed/molden.inp'], seed['files']
entry = seed['files'][0]
assert set(entry) == {'path', 'bytes', 'sha256'}, entry
assert entry['bytes'] == (dalton / 'molden.inp').stat().st_size
assert entry['sha256'] == sha256(dalton / 'molden.inp')
assert seed['source_files'][0]['source'] == str(dalton / 'molden.inp'), seed['source_files']
assert seed['method'] == 'unknown' and 'could not be checked' in seed['method_check'], seed
seeded_from = seed['seeded_from']
for key in ('seed', 'stage', 'molden', 'dalton_dir', 'n_occ', 'n_ao', 'n_mo', 'L', 'thresh',
            'max_offdiag_pre_loewdin', 'archive', 'preserved_copy'):
    assert key in seeded_from, (key, seeded_from)
assert seeded_from['n_occ'] == 5 and seeded_from['stage'] == 'ground_state', seeded_from
copy = work / 'Seeded/gs/checkpoints' / generation / 'seed/molden.inp'
assert sha256(copy) == entry['sha256']
settings = m['metrics']['convergence_report']['convergence_settings']
assert settings['seed'] == 'Dalton directory ' + str(dalton) + ' (user)', settings.get('seed')
cold_generation, cold_manifest = manifest(work / 'Cold')
assert not cold_manifest.get('seed'), cold_manifest.get('seed')
assert 'seed' not in cold_manifest['metrics']['convergence_report']['convergence_settings']

# A saved state wins: --continue loads it and seeds nothing.
again = work / 'SeededAgain'
shutil.copytree(work / 'Seeded', again)
log = run('seeded-continue', again, seeded, resume=True)
record = result(again)
assert record['work'] == 'loaded' and record['iterations_this_run'] == 0, (record['work'], record['iterations_this_run'])
assert '[DALTON-SEED]' not in log, [line for line in log.splitlines() if 'DALTON' in line]

# An altered seed copy is refused when the generation is read.
tampered = work / 'Tampered'
shutil.copytree(work / 'Seeded', tampered)
target = tampered / 'gs/checkpoints' / generation / 'seed/molden.inp'
target.chmod(0o644)
target.write_text(target.read_text() + '\n')
log = run('tampered-continue', tampered, seeded, resume=True, success=False)
assert 'Dalton seed files: seed/molden.inp altered' in log, log[-3000:]

# Another geometry is refused by the seed's fingerprint before the SCF.
moved = seeded.replace('0.1925923500', '0.2925923500')
assert moved != seeded
log = run('moved', work / 'Moved', moved, success=False)
assert 'GEOMETRY FINGERPRINT MISMATCH' in log, log[-3000:]
assert not (work / 'Moved/gs/checkpoints/latest').exists()

print(json.dumps({'cold_iterations': cold_iterations, 'seeded_iterations': seeded_iterations,
                  'energy_cold': energy_c, 'energy_seeded': energy_s,
                  'mpi_ranks': a.ranks}, indent=2))
print('PASS')
