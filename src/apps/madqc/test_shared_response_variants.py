#!/usr/bin/env python3
"""Native shared records for existing open-shell driven solver variants.

Loose gates isolate wiring and payload shape. Scientific accuracy is checked
separately at the end of response integration, with unchanged references.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--atom', choices=['H', 'F', 'C'], default='H')
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--nuclear', action='store_true')
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
work = Path(tempfile.mkdtemp(prefix='response-variants-'))
print('artifacts:', work, flush=True)
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
binary = str(Path(a.binary).resolve())
na, nb, nopen = (1, 0, 0) if a.nuclear else {'H': (1, 0, 1), 'F': (5, 4, 1), 'C': (4, 2, 2)}[a.atom]
atom = 'He' if a.nuclear else a.atom
nuclear = 'nuclear true\n nuclear.directions z\n nuclear.frequencies [0.0]' if a.nuclear else ''
subjects = [('dipole_z', 0.), ('dipole_z', .05)] + ([('nuc_0_z', 0.)] if a.nuclear else [])
(work/'input.in').write_text(f'''dft
 xc hf
 localize canon
 spin_restricted {str(a.nuclear).lower()}
 nopen {nopen}
 protocol [1e-4]
 l 30
 maxiter 1
 convergence.bsh_residual 1000
 convergence.density_change 1000
 convergence.energy_change 1000
 print_level 0
end
response
 dipole true
 dipole.directions z
 dipole.frequencies [0.0, 0.05]
 requested_properties [polarizability]
 {nuclear}
 convergence.bsh_residual 1000
 convergence.density_change 1000
 maxiter 1
 print_level 1
end
molecule
 units atomic
 no_orient true
 {atom} 0 0 0
end
''')

def run(label, extra=(), success=True):
    result = subprocess.run(launch+[binary, 'input.in', '--wf=response', '--prefix=run']+list(extra),
        cwd=work, env=env, capture_output=True, text=True, timeout=240)
    (work/f'{label}.log').write_text(result.stdout+result.stderr)
    assert (result.returncode == 0) == success, (result.returncode, result.stdout[-5000:], result.stderr)

root = work/'run'
def records():
    result = []
    for perturbation, frequency in subjects:
        task = root/'checkpoints/response/tasks'/f'fd:{perturbation}@f{frequency:.5f}'
        generation = json.loads((task/'latest').read_text())['generation']
        result.append(json.loads((task/generation/'manifest.json').read_text()))
    return result

run('initial')
first = records()
for (_, frequency), record in zip(subjects, first):
    assert record['task']['outcome']['complete']
    metrics = {m['metric']: m for m in record['task']['iteration']['metrics']}
    assert metrics['energy_change']['status'] == 'not_applicable'
    # One update with loose targets: change metrics have a predecessor.
    assert record['task']['iteration']['updates_this_attempt'] == 1
    assert metrics['density_change']['status'] == 'measured'
    assert metrics['kinetic_change']['status'] == 'measured'
    expected = {f'x_alpha:{i}' for i in range(na)} | {f'x_beta:{i}' for i in range(nb)}
    if frequency:
        expected |= {f'y_alpha:{i}' for i in range(na)} | {f'y_beta:{i}' for i in range(nb)}
    assert {c['id'] for c in metrics['bsh_residual']['evidence']['components']} == expected
response = next(c for c in json.loads((root/'results.json').read_text())['calculations']
                if c['calculation_id'] == 'response')
for tensor in response['properties']['polarizability']['values']:
    if a.nuclear:
        assert tensor['status'] == 'complete'
    else:
        assert tensor['status'] == 'unavailable' and not tensor['components']
        assert 'not implemented' in tensor['availability']['zz']
run('reuse')
assert records() == first
run('continue', ['--response=convergence.bsh_residual=0;maxiter=1;kain=true'], success=False)
for record in records():
    assert record['task']['iteration']['updates_this_attempt'] == 1
    metrics = {m['metric']: m for m in record['task']['iteration']['metrics']}
    assert metrics['density_change']['status'] == 'measured'
    assert metrics['kinetic_change']['status'] == 'measured'
print(f'{atom}: static/dynamic shared records, reuse and update passed; nuclear={a.nuclear}', flush=True)
