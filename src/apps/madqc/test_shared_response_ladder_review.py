#!/usr/bin/env python3
"""Regression checks for frequency identity, archive layout and current result mirrors."""
from response_test_checks import assert_same_alpha_properties
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser();p.add_argument('--binary',required=True)
p.add_argument('--ranks',type=int,default=1);p.add_argument('--mpiexec',default='mpiexec')
a=p.parse_args();work=Path(tempfile.mkdtemp(prefix='response-ladder-review-'));print('artifacts:',work,flush=True)
(work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 maxiter 25
 convergence.bsh_residual 1e-4
 convergence.density_change 1e-4
 convergence.energy_change 1e-6
 print_level 0
end
response
 dipole.directions z
 dipole.frequencies [0.05,0.1]
 convergence.bsh_residual 1000
 convergence.density_change 1000
 maxiter 1
 print_level 1
end
molecule
 units atomic
 He 0 0 0
end
''')
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
launch=[] if a.ranks==1 else [a.mpiexec,'-n',str(a.ranks)]
command=launch+[str(Path(a.binary).resolve()),'input.in','--wf=response','--prefix=run']
root=work/'run';tasks=root/'checkpoints/response/tasks'
metadata=work/'run/analysis/response/response_metadata.json'
failures=[]
def run(label,options=''):
    r=subprocess.run(command+(['--response='+options] if options else []),cwd=work,env=env,capture_output=True,text=True,timeout=180)
    (work/(label+'.log')).write_text(r.stdout+r.stderr)
    return r.returncode

def saved(freq):
    t=tasks/f'fd:dipole_z@f{freq:.5f}';g=json.loads((t/'latest').read_text())['generation']
    return g,json.loads((t/g/'manifest.json').read_text())
def result():return next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
def mirror_frequencies():return {r['omega'] for rows in json.loads(metadata.read_text())['properties']['alpha'].values() for r in rows}
def verify_reuse(before,after,original):
    assert after['attempt_id']!=before['attempt_id'],'result snapshot was not republished'
    assert_same_alpha_properties(before['properties'],after['properties'])
    sources=json.loads((root/'checkpoints/response/attempts'/after['attempt_id']/'result_sources.json').read_text())['states']
    assert set(sources)==set(original)
    for task,(generation,manifest) in original.items():
        assert sources[task]['generation']==generation and sources[task]['attempt_id']==manifest['attempt_id']

assert run('initial')==0
original={f'fd:dipole_z@f{f:.5f}':saved(f) for f in [.05,.1]}
before=result();metadata.unlink()
assert run('reuse')==0
verify_reuse(before,result(),original)
try:
    verify_reuse(before,before,original)
except AssertionError:
    pass
else:raise AssertionError('reuse guard accepted a stale results snapshot')
pointer=(tasks/'fd:dipole_z@f0.05000/latest').read_bytes()
payloads={f:hashlib.sha256(f.read_bytes()).hexdigest() for f in tasks.glob('*/*/*.00000')}
code=run('alias','dipole.frequencies=[0.050001]')
if code==0 or (tasks/'fd:dipole_z@f0.05000/latest').read_bytes()!=pointer or any(not f.exists() or hashlib.sha256(f.read_bytes()).hexdigest()!=h for f,h in payloads.items()):
    failures.append('cross-invocation frequency alias accepted or altered saved state')
analysis=json.loads(metadata.read_text())
analysis['properties']['unrelated_analysis']={'marker':'preserve'}
metadata.write_text(json.dumps(analysis))
assert run('loose-again')==0
residuals={f:next(m['value'] for m in saved(f)[1]['task']['iteration']['metrics'] if m['metric']=='bsh_residual') for f in [.05,.1]}
gate=sum(residuals.values())/2
assert min(residuals.values())<gate<max(residuals.values())
# No update: the deck allows one, which could carry the unmet state under the gate.
assert run('partial',f'convergence.bsh_residual={gate:.17g};maxiter=0')!=0
wanted={f for f,v in residuals.items() if v<=gate}
if mirror_frequencies()!=wanted:failures.append('analysis mirror retained an unqualified alpha tensor')
assert run('subset','dipole.frequencies=[0.1]')==0
if mirror_frequencies()!={.1}:failures.append('analysis mirror retained an unrequested frequency')
assert json.loads(metadata.read_text())['properties']['unrelated_analysis']=={'marker':'preserve'}
info=json.loads((work/'run.calc_info.json').read_text())
properties=info['tasks'][1]['properties']['response_properties']
if set(properties)!={'alpha'}:failures.append('calc_info exported unrelated analysis as requested properties')
exported=properties['alpha']
if {r['omega'] for rows in exported.values() for r in rows}!={.1}:failures.append('calc_info exported an unrequested frequency')
if run('small-frequency','dipole.frequencies=[1e-13]')!=0:
    failures.append('small nonzero static frequency could not load its native archive')
else:
    assert saved(1e-13)[1]['task']['definition']['settings']['layout']=='closed_shell_static'
assert not failures,failures
print('frequency identity/layout, mirror and reuse-source checks passed')
