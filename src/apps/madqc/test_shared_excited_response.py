#!/usr/bin/env python3
"""Coupled bundle lifecycle; loose gates test wiring, not excitation accuracy."""
from response_test_checks import assert_energy_only_es_properties
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--full',action='store_true')
p.add_argument('--ranks',type=int,default=1)
p.add_argument('--mpiexec',default='mpiexec')
a=p.parse_args()
work=Path(tempfile.mkdtemp(prefix='shared-excited-'))
print('artifacts:',work,flush=True)
(work/'input.in').write_text(f'''dft
 xc hf
 localize canon
 protocol [1e-4,1e-5]
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
 excited.enable true
 excited.tda {str(not a.full).lower()}
 excited.num_states 2
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 1
end
molecule
 units angstrom
 no_orient true
 H 0 0 -0.37
 H 0 0 0.37
end
''')
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
cmd=([] if a.ranks==1 else [a.mpiexec,'-n',str(a.ranks)])+[str(Path(a.binary).resolve()),'input.in','--wf=response','--prefix=run']
root=work/'run'
model='full' if a.full else 'tda'
task=root/'checkpoints/response/tasks'/f'es:{model}'
def run(label,options='',success=True):
 r=subprocess.run(cmd+(['--response='+options] if options else []),cwd=work,env=env,capture_output=True,text=True,timeout=300)
 (work/f'{label}.log').write_text(r.stdout+r.stderr)
 assert (r.returncode==0)==success,(label,r.returncode,r.stdout[-4500:],r.stderr)
def saved():
 generation=json.loads((task/'latest').read_text())['generation']
 return json.loads((task/generation/'manifest.json').read_text())
run('initial')
first=saved()
assert first['task']['definition']['protocol']=={'threshold':1e-5,'k':8}
assert first['task']['outcome']['complete']
assert len(first['archives'])==2
assert set(first['task']['definition']['subjects'])=={'es_root_0000','es_root_0001'}
for m in first['task']['iteration']['metrics']:
 if m['metric'] in ['energy_change','density_change','kinetic_change']:
  assert m['status']=='unavailable'
run('reuse')
assert saved()==first
run('tighten','convergence.bsh_residual=0;maxiter=1',success=False)
tight=saved()
assert tight['task']['iteration']['updates_this_attempt']==1
assert not tight['task']['outcome']['complete']
run('add-root','excited.num_states=3')
expanded=saved()
assert len(expanded['archives'])==3
assert expanded['task']['definition']['protocol']=={'threshold':1e-5,'k':8}
assert expanded['task']['iteration']['updates_this_attempt']==0
for m in expanded['task']['iteration']['metrics']:
 if m['metric']=='energy_change': assert m['status']=='unavailable'
run('expanded-reuse','excited.num_states=3')
assert saved()==expanded
result=next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
assert 'excitation_energy' in result['properties']
env['MADQC_ANALYZE_PROBE']='0'
run('diagnostic-probe','excited.num_states=3')
assert saved()==expanded
result=next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
assert_energy_only_es_properties(result['properties'])
print('shared coupled bundle restart checks passed',flush=True)
