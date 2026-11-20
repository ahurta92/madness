#!/usr/bin/env python3
"""FD, ES and VBC native inspection must leave all saved run files unchanged."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--solver')
p.add_argument('--run-root',action='append',default=[])
a=p.parse_args()
work=Path(tempfile.mkdtemp(prefix='shared-inspection-'))
print('artifacts:',work,flush=True)
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
roots=[Path(r).resolve() for r in a.run_root]
if not roots:
 assert a.solver
 (work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 l 30
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 0
end
response
 dipole true
 dipole.directions z
 dipole.frequencies [0.0]
 quadratic true
 excited.enable true
 excited.tda true
 excited.num_states 1
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 0
end
molecule
 units angstrom
 no_orient true
 H 0 0 -0.37
 H 0 0 0.37
end
''')
 r=subprocess.run([str(Path(a.solver).resolve()),'input.in','--wf=response','--prefix=run'],cwd=work,env=env,capture_output=True,text=True,timeout=300)
 (work/'solve.log').write_text(r.stdout+r.stderr)
 assert r.returncode==0,(r.stdout[-2500:],r.stderr)
 roots=[work/'run']
def hashes(root):
 return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file()}
remaining={'fd:','es:','vbc:'}
for root in roots:
 before=hashes(root)
 for task in (root/'checkpoints/response/tasks').iterdir():
  family=next((k for k in remaining if task.name.startswith(k)),None)
  if not family: continue
  generation=json.loads((task/'latest').read_text())['generation']
  if generation is None: continue
  output=work/family[:-1]
  r=subprocess.run([str(Path(a.binary).resolve()),'--checkpoint-root='+str(root),'--task='+task.name,'--out='+str(output),'--max-orbitals=1'],cwd=work,env=env,capture_output=True,text=True,timeout=120)
  (work/(family[:-1]+'.log')).write_text(r.stdout+r.stderr)
  assert r.returncode==0,(task.name,r.stdout[-2500:],r.stderr)
  report=json.loads((output/'inspection.json').read_text())
  assert report['source']['generation']==generation
  assert list(output.glob('*_recon.json'))
  assert list(output.glob('*_comp.json'))
  if family!='vbc:': assert report['metrics']
  remaining.remove(family)
 assert hashes(root)==before,'inspection mutated the saved calculation'
assert not remaining,remaining
print('native inspection checks passed',flush=True)
