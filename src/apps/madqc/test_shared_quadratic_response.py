#!/usr/bin/env python3
"""Native VBC dependency/reuse wiring, without a numerical beta reference claim."""
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
p.add_argument('--subworlds', type=int, default=0)
a = p.parse_args()
work = Path(tempfile.mkdtemp(prefix='shared-quadratic-'))
print('artifacts:', work, flush=True)
(work/'input.in').write_text('''dft
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
 dipole.frequencies [0.02]
 quadratic true
 beta.frequencies [0.02]
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 1
end
molecule
 units angstrom
 no_orient true
 He 0 0 0
end
''')
if a.subworlds:
    deck=work/'input.in'
    deck.write_text(deck.read_text().replace(' dipole.directions z',' dipole.directions xz').replace('response\n','response\n subworlds '+str(a.subworlds)+'\n'))
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', OMPI_MCA_hwloc_base_binding_policy='none')
cmd = ([] if a.ranks==1 else [a.mpiexec,'-n',str(a.ranks)]) + [str(Path(a.binary).resolve()), 'input.in', '--wf=response', '--prefix=run']
root = work/'run'
tasks = root/'checkpoints/response/tasks'
def run(label, options='', success=True):
    r = subprocess.run(cmd+(['--response='+options] if options else []), cwd=work, env=env, capture_output=True, text=True, timeout=300)
    (work/f'{label}.log').write_text(r.stdout+r.stderr)
    assert (r.returncode==0)==success, (label,r.returncode,r.stdout[-4000:],r.stderr)
def manifest(path):
    gen = json.loads((path/'latest').read_text())['generation']
    return json.loads((path/gen/'manifest.json').read_text())
def response():
    return next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
def beta():
    return response()['properties']['hyperpolarizability']['values'][0]
run('initial')
legacy=json.loads((work/'run.calc_info.json').read_text())['tasks'][1]['properties']['response_properties']
assert legacy['beta']['1e-05_k8'], legacy
vbc = next(p for p in tasks.iterdir() if p.name.startswith('vbc:'))
first = manifest(vbc)
assert first['task']['kind']=='direct'
assert 'iteration' not in first['task']
assert first['task']['definition']['protocol']=={'threshold':1e-5,'k':8}
assert {x['role'] for x in first['task']['inputs']}=={'ground','B','C'}
assert len(first['archives'])==1
assert beta()['status']=='complete' and 'zzz' in beta()['components']
run('reuse')
assert manifest(vbc)==first
run('tighten','convergence.bsh_residual=0;maxiter=1',False)
assert manifest(vbc)==first # invalidated inputs never certify old VBC or properties
assert beta()['status']=='unavailable' and not beta()['components']
run('requalify')
second = manifest(vbc)
assert second['generation']!=first['generation']
assert second['task']['inputs']!=first['task']['inputs']
assert beta()['status']=='complete'
run('reuse-rebuilt')
assert manifest(vbc)==second
attempt = json.loads((root/'timings.json').read_text())['attempts'][-1]
def descendants(node):
    yield node
    for child in node.get('children',[]): yield from descendants(child)
nodes = list(descendants(attempt['timing']))
assert any(n['name']=='vbc' for n in nodes)
assert not any(n['name']=='build_source' for n in nodes)
assert any(n['name']=='contraction' for n in nodes)
if a.subworlds:
    from response_test_checks import assert_subworld_timing_scopes
    assert_subworld_timing_scopes(json.loads((root/'timings.json').read_text())['attempts'], a.subworlds, a.ranks)
print('shared quadratic dependency checks passed', flush=True)
