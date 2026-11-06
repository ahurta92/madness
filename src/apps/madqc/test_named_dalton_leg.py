#!/usr/bin/env python3
"""Driven legs seeded from a Dalton directory on the named path (D-131): water
at the geometry of the qctest response_h2o_tpa_seeded, alpha_zz at that
case's two DIPLEN frequencies (its 2PA legs, omega_f/2, as the RSPVEC lists
them: the import seeds a leg to disk only on an exact match), seeded from its
dalton/. The import writes each leg's seed under the step
"<step>_dseed"; the stage saves with them, every generation holds the seed
archives (the generation published after the last leg too, so a --continue
interrupted between legs could still start from them), the seeds are neither
state records nor report channels, and --continue twice reuses every leg. A
seed entry marked converged is refused on read."""
import argparse
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
work = named_workdir.work_directory('named-dalton-leg-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
dalton = (Path(a.source).parent.parent / 'examples/qc/response_h2o_tpa_seeded/dalton').resolve()
assert (dalton / 'RSPVEC').exists(), dalton

deck = """molecule water
  eprec 1e-6
  units atomic
  no_orient true
  O   0.0000000000   0.0000005943   0.1925923500
  H   0.0000000000   1.4214240411  -0.8721230117
  H   0.0000000000  -1.4214246355  -0.8721234678
end

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

molresponse alpha
  reference gs
  protocol [1e-4]
  dconv 1e-4
  maxiter 30
  requested_properties [polarizability]
  polarizability.frequencies [0.160930348, 0.191670362]
  polarizability.frequency_units au
  polarizability.components [zz]
  seed.dalton_dir """ + str(dalton) + """
  seed.freq_tol 0.01
end
"""


def run(label, root, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(deck)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=1500)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


def generation(root):
    store = root / 'alpha/checkpoints'
    name = json.loads((store / 'latest').read_text())['generation']
    return store / name, json.loads((store / name / 'manifest.json').read_text())


root = work / 'Seeded'
log = run('seeded', root)
stage = named_records.stage_record(root / 'alpha')
assert stage['status'] == 'completed', stage['status']
# Both legs started from their Dalton seed (the engine's loader reads the seed
# step as an inexact source and refines it).
starts = [line for line in log.splitlines()
          if line.startswith('[LOAD] try_load_fd_state') and '_dseed' in line and 'exact= false' in line]
assert len(starts) == 2, [line for line in log.splitlines() if 'LOAD' in line or 'DALTON' in line]

# The generation published after the last leg still holds both seeds, and the
# cleanup of older generations left them on disk in it.
directory, m = generation(root)
metadata = m['numerical']['engine_metadata']
seed_steps = [k for k in metadata['protocols'] if k.endswith('_dseed')]
assert len(seed_steps) == 1, metadata['protocols']
step = seed_steps[0][:-len('_dseed')]
seeds = metadata['fd_states']['dipole_z'][seed_steps[0]]
assert len(seeds) == 2 and all(e['seed'] == 'dalton_import' and not e['converged'] and e['iter'] == 0
                               for e in seeds.values()), seeds
solved = metadata['fd_states']['dipole_z'][step]
assert len(solved) == 2 and all(e['converged'] for e in solved.values()), solved
archives = {x['basename'] for x in m['archives']}
for e in seeds.values():
    assert e['archive'] in archives, (e['archive'], sorted(archives))
    assert any(directory.rglob(e['archive'] + '.*')), (e['archive'], directory)
store = root / 'alpha/checkpoints'
print('generations on disk:', sorted(p.name for p in store.iterdir() if p.is_dir()), flush=True)

# Seeds are not solves: no state record and no report channel.
assert all('_dseed' not in s['key'] for s in m['numerical']['states']), [s['key'] for s in m['numerical']['states']]
channels = stage['convergence_report']['channels']
assert len(channels) == 2 and all(c['protocol_step']['thresh'] == metadata['protocols'][step]['thresh']
                                  and '_dseed' not in c['archive'] for c in channels), channels
alpha = stage['properties']['polarizability']

# --continue from that generation, twice: every leg is reused.
for label in ('continue-1', 'continue-2'):
    run(label, root, resume=True)
    again = named_records.stage_record(root / 'alpha')
    assert again['work'] == 'loaded' and again['iterations_this_run'] == 0, (again['work'], again['iterations_this_run'])
    assert again['properties']['polarizability'] == alpha

# A seed entry marked converged is refused when the generation is read.
tampered = work / 'Tampered'
shutil.copytree(root, tampered)
target, tm = generation(tampered)
point = next(iter(tm['numerical']['engine_metadata']['fd_states']['dipole_z'][seed_steps[0]]))
tm['numerical']['engine_metadata']['fd_states']['dipole_z'][seed_steps[0]][point]['converged'] = True
manifest = target / 'manifest.json'
manifest.chmod(0o644)
manifest.write_text(json.dumps(tm, indent=2))
log = run('tampered-continue', tampered, resume=True, success=False)
assert 'invalid response state: Dalton seed' in log, log[-3000:]

print(json.dumps({'seed_step': seed_steps[0], 'iterations': stage['iterations_this_run'],
                  'mpi_ranks': a.ranks}, indent=2))
print('PASS')
