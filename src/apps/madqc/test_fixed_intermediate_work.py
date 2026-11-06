#!/usr/bin/env python3
"""Fixed intermediate work reaches the final representation in named solvers."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records
import named_workdir

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--source', required=True)
a = parser.parse_args()
binary = str(Path(a.binary).resolve())
source = Path(a.source)
work = named_workdir.work_directory('fixed-intermediate-work-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '4')

# Loose checks deliberately pass before the eight-iteration coarse allowance.
# Use the maintained named H2 fixture; these are control-flow assertions,
# not reference values for energies or response properties.
base = (source / 'named_h2.in').read_text()
head, optical = base.split('molresponse optical', 1)
base = head.replace('  maxiter 30\n', '  maxiter [8, 30]\n') + 'molresponse optical' + optical.replace(
    '  maxiter 30\n', '  maxiter 8\n')
base = base.replace('  dconv 1e-5\n', '  dconv 1e-3\n')
base = base.replace('  dconv 1e-4\n', '  dconv 1e-1\n')
_, es = (source / 'named_h2_es_tda.in').read_text().split('molresponse excited', 1)
es = 'molresponse excited' + es
es = es.replace('  protocol [1e-4]\n', '  protocol [1e-4, 1e-6]\n')
es = es.replace('  maxiter 20\n', '  maxiter 8\n')
es = es.replace('  dconv 1e-4\n', '  dconv 1e-1\n')
es = es.replace('[excited_states, polarizability]', '[excited_states]')
es = '\n'.join(line for line in es.split('\n') if 'polarizability.' not in line)
deck = base + '\n' + es + '\n' + es.replace('molresponse excited', 'molresponse full').replace(
    'excited_states.method tda', 'excited_states.method rpa')
inp = work / 'input.in'
inp.write_text(deck)
root = work / 'run'

def run(label, resume=False):
    log = work / (label + '.log')
    with log.open('w') as out:
        result = subprocess.run(
            [binary, str(inp), '--prefix=' + str(root)] + (['--continue'] if resume else []),
            cwd=work, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=900)
    assert result.returncode == 0, log.read_text()[-8000:]
    return named_records.stages(root)

def manifest(stage):
    store = root / stage / 'checkpoints'
    generation = json.loads((store / 'latest').read_text())['generation']
    return json.loads((store / generation / 'manifest.json').read_text())

stages = run('fresh')
assert len(stages) == 4
history = stages[0]['convergence_report']['iteration_history']
coarse = [row for row in history if row['kinetic_trace']['sample']['threshold'] == 1e-4]
fine = [row for row in history if row['kinetic_trace']['sample']['threshold'] == 1e-6]
assert len(coarse) == 8 and fine, (len(coarse), len(fine))
assert any(row['configured_progress']['fixed_intermediate_work'] for row in coarse)
for name in ('optical', 'excited', 'full'):
    records = manifest(name)['numerical']['states']
    coarse = [record for record in records if record['protocol_step']['thresh'] == 1e-4]
    fine = [record for record in records if record['protocol_step']['thresh'] == 1e-6]
    assert coarse and fine, (name, records)
    assert all(record['iterations_at_protocol_step'] == 8 for record in coarse), (name, coarse)
    assert all(record['iterations_at_protocol_step'] > 0 for record in fine), (name, fine)
    if name != 'optical':
        # The loose one-root case locks before the fixed coarse count ends.
        # Later all-locked steps must retain comparable eigenvalue evidence.
        assert all(record['stop_reason'] == 'converged' for record in coarse), (name, coarse)
        for record in coarse:
            assert all(math.isfinite(root['excitation_energy_change'])
                       for root in record['roots']), (name, record)
        metadata = manifest(name)['numerical']['engine_metadata']
        locked_iterations = {
            sample['kinetic_energy']['sample']['iteration']
            for entry in metadata['excited_states'].values()
            for sample in entry['kinetic_history']
            if sample['kinetic_energy']['sample']['threshold'] == 1e-4
            and sample['root_locked']
        }
        assert len(locked_iterations) >= 2, (name, locked_iterations)
    print(name, 'fixed coarse work and final protocol verified', flush=True)
assert all(stage['status'] == 'completed' for stage in stages)
continued = run('continued', resume=True)
assert all(stage['iterations_this_run'] == 0 for stage in continued), continued
print('fixed intermediate work acceptance passed', flush=True)
