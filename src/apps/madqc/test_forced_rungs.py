#!/usr/bin/env python3
"""Explicit diagnostic counts must not stop early at ordinary convergence."""
import argparse
import json
import sys
import os
from pathlib import Path
import subprocess
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
work = named_workdir.work_directory('forced-rungs-')
print('artifacts:', work, flush=True)
binary = str(Path(a.binary).resolve())
deck = (Path(a.source) / 'named_h2.in').read_text()
deck = deck.replace('moldft gs\n', 'moldft gs\n  diagnostic.iterations_per_rung 5\n')
deck = deck.replace('molresponse optical\n', 'molresponse optical\n  diagnostic.iterations_per_rung 8\n')
input_file = work / 'forced.in'
input_file.write_text(deck)
root = work / 'run'
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
with (work / 'run.log').open('w') as out:
    result = subprocess.run(launch + [binary, str(input_file), '--prefix=' + str(root)],
                            cwd=work, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=600)
assert result.returncode == 0, (work / 'run.log').read_text()[-5000:]
stages = named_records.stages(root)
assert stages[0]['iterations_this_run'] == 10, stages[0]['iterations_this_run']
assert stages[1]['iterations_this_run'] == 32, stages[1]['iterations_this_run']
for channel in stages[1]['convergence_report']['channels']:
    assert channel['iterations_at_save'] == 8
    assert len(channel['density_history']) == 8
    assert len(channel['function_kinetic_history']) == 9
assert len(stages[0]['convergence_report']['iteration_history']) == 10
print('forced per-rung diagnostics passed', flush=True)
