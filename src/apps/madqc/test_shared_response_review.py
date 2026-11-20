#!/usr/bin/env python3
"""Rejected alpha requests must not prepare response or change saved provenance."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
work = Path(tempfile.mkdtemp(prefix='response-review-'))
print('artifacts:', work, flush=True)
deck = '''dft
 xc hf
 localize canon
 k 6
 protocol [1e-4]
 maxiter 25
 convergence.criteria [bsh_residual, density_change, energy_change]
 convergence.bsh_residual 1e-4
 convergence.density_change 1e-4
 convergence.energy_change 1e-6
 print_level 0
end
response
 dipole.directions z
 dipole.frequencies [0.05]
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 maxiter 0
 print_level 1
end
molecule
 units atomic
 He 0 0 0
end
'''
(work/'input.in').write_text(deck)
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]


def run(label, extra=()):
    r = subprocess.run(launch+[str(Path(a.binary).resolve()), 'input.in', '--wf=response',
                       '--prefix=run']+list(extra), cwd=work, env=env,
                       capture_output=True, text=True, timeout=180)
    text = r.stdout+r.stderr
    (work/(label+'.log')).write_text(text)
    return r.returncode, text


assert run('initial')[0] == 0
store = work/'run/checkpoints'
task = 'fd:dipole_z@f0.05000'
latest = store/'response/tasks'/task/'latest'
original = latest.read_bytes()
failures = []
for label, value in [('property', 'requested_properties=[polarizability,bogus]'),
                     ('axis', 'dipole.directions=xq'), ('subworld', 'subworlds=1'),
                     ('frequency-alias', 'dipole.frequencies=[0.05,0.050001]')]:
    code, text = run(label, ['--response='+value])
    if code == 0 or 'GROUND_STATE_INFO' in text or latest.read_bytes() != original:
        failures.append(label+': accepted or prepared an unsupported response request')


def inputs():
    calculation = next(c for c in json.loads((store/'run.json').read_text())['calculations']
                       if c['calculation_id'] == 'response')
    return next(t['inputs'] for t in calculation['tasks'] if t['task_id'] == task)


before = inputs()
code, text = run('changed-ground', ['--dft=localize=boys'])
assert code != 0 and ('model' in text or 'source' in text or 'inputs' in text), text[-4000:]
if inputs() != before or latest.read_bytes() != original:
    failures.append('changed-source rejection altered registered provenance')
assert not failures, failures
print('shared response rejection and provenance checks passed')
