#!/usr/bin/env python3
"""Shared Raman/2PA lifecycle; deliberately loose gates isolate the integration."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--raman',action='store_true')
p.add_argument('--full-raman',action='store_true')
a=p.parse_args()
if a.full_raman: a.raman=True
work=Path(tempfile.mkdtemp(prefix='shared-derived-'))
print('artifacts:',work,flush=True)
request=''' requested_properties [raman]
 raman.atoms [0]
 raman.directions z
 raman.frequencies [0.0]
''' if a.raman else ''' excited.enable true
 excited.tda true
 excited.num_states 1
 excited.tpa true
 convergence.energy_change 1000
'''
(work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4,1e-5]
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
 dipole.frequencies [0.0]
 maxiter 1
 convergence.bsh_residual 1000
 convergence.density_change 1000
 print_level 0
'''+request+'''end
molecule
 units angstrom
 no_orient true
 H 0 0 -0.37
 H 0 0 0.37
end
''')
if a.full_raman:
 deck=work/'input.in'
 deck.write_text(deck.read_text().replace(' raman.atoms [0]\n raman.directions z',' raman.nuc_atom -1').replace(' dipole.directions z',' dipole.directions xyz'))
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
cmd=[str(Path(a.binary).resolve()),'input.in','--wf=response','--prefix=run']
root=work/'run'
def run(label,options='',success=True):
 r=subprocess.run(cmd+(['--response='+options] if options else []),cwd=work,env=env,capture_output=True,text=True,timeout=480)
 (work/f'{label}.log').write_text(r.stdout+r.stderr)
 assert (r.returncode==0)==success,(label,r.returncode,r.stdout[-3000:],r.stderr)
def result():
 return next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
def manifests():
 records={}
 for path in (root/'checkpoints/response/tasks').iterdir():
  if not path.is_dir() or not (path/'latest').is_file(): continue
  generation=json.loads((path/'latest').read_text())['generation']
  if generation and (path/generation/'manifest.json').is_file():
   records[path.name]=json.loads((path/generation/'manifest.json').read_text())
 return records
run('initial')
legacy=json.loads((work/'run.calc_info.json').read_text())['tasks'][1]['properties']['response_properties']
assert legacy['raman' if a.raman else 'tpa'], legacy
first=manifests()
prop='polarizability_gradient' if a.raman else 'two_photon_absorption'
assert result()['properties'][prop]['values']
assert all(row['status']=='complete' for row in result()['properties'][prop]['values'])
if a.raman and not a.full_raman:
 assert 'vibrational_analysis' not in result()['properties'] # partial coordinate request
run('reuse')
assert manifests()==first
if a.full_raman:
 properties=result()['properties']
 assert len(properties['vibrational_analysis']['values']['hessian'])==6
 assert not properties['vibrational_analysis']['availability']
 assert not properties['raman_spectrum']['availability']
 print('shared full Raman native reuse checks passed',flush=True)
 raise SystemExit(0)
run('tighten','convergence.bsh_residual=0;maxiter=1',False)
assert all(row['status']=='unavailable' for row in result()['properties'][prop]['values'])
run('requalify')
second=manifests()
assert all(row['status']=='complete' for row in result()['properties'][prop]['values'])
if not a.raman:
 derived=[m for m in second.values() if 'derived_from' in m['task']['definition'].get('settings',{})]
 assert len(derived)>=3
 for m in derived:
  assert m['task']['definition']['protocol']=={'threshold':1e-5,'k':8}
  # New states at the requalified root frequency: one update meets the loose targets.
  assert m['task']['iteration']['updates_this_attempt']==1
run('reuse-refined')
assert manifests()==second
print('shared derived response checks passed',flush=True)
