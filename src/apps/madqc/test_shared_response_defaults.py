#!/usr/bin/env python3
"""Default response selects the shared contract; calibrated values are checked separately."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--families',action='store_true')
a=p.parse_args()
work=Path(tempfile.mkdtemp(prefix='shared-defaults-'))
print('artifacts:',work,flush=True)
(work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 l 30
 maxiter 30
 convergence.bsh_residual 1e-4
 convergence.density_change 1e-4
 convergence.energy_change 1e-6
 print_level 0
end
response
 dipole true
 dipole.directions z
 dipole.frequencies [0.0]
 maxiter 30
 print_level 1
end
molecule
 units atomic
 He 0 0 0
end
''')
if a.families:
 (work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 l 30
 maxiter 1
 convergence.bsh_residual 1000
 convergence.density_change 1000
 convergence.energy_change 1000
 print_level 0
end
response
 requested_properties [polarizability, hyperpolarizability, raman]
 dipole.directions z
 dipole.frequencies [0.0]
 beta.frequencies [0.0]
 raman.atoms [0]
 raman.directions z
 raman.frequencies [0.0]
 excited.enable true
 excited.tda true
 excited.num_states 1
 excited.tpa true
 maxiter 3
 convergence.bsh_residual 1000
 convergence.density_change 1000
 convergence.energy_change 1000
 print_level 0
end
molecule
 units angstrom
 no_orient true
 H 0 0 -0.37
 H 0 0 0.37
end
''')
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
r=subprocess.run([str(Path(a.binary).resolve()),'input.in','--wf=response','--prefix=run'],cwd=work,env=env,capture_output=True,text=True,timeout=600)
(work/'initial.log').write_text(r.stdout+r.stderr)
assert r.returncode==0,(r.stdout[-2000:],r.stderr)
root=work/'run'
response=next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
assert response['properties']['polarizability']['values'][0]['status']=='complete'
task=root/'checkpoints/response/tasks/fd:dipole_z@f0.00000'
generation=json.loads((task/'latest').read_text())['generation']
record=json.loads((task/generation/'manifest.json').read_text())
assert record['task']['outcome']['complete']
assert record['task']['iteration']['updates_this_attempt']>0
assert (root/'timings.json').is_file()
if a.families:
 for name in ['polarizability','hyperpolarizability','polarizability_gradient','two_photon_absorption']:
  assert response['properties'][name]['values'],name
  assert all(row['status']=='complete' for row in response['properties'][name]['values']),name
 assert response['properties']['excitation_energy']['components']
 tasks=list((root/'checkpoints/response/tasks').iterdir())
 for family in ['fd:','es:','vbc:']:
  assert any(task.name.startswith(family) for task in tasks),family
print('shared response default routing passed',flush=True)
