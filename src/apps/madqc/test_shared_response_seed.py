#!/usr/bin/env python3
"""DALTON functions enter native shared ES checkpoints as freshly measured guesses."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--source',required=True)
a=p.parse_args()
case=Path(a.source).resolve()/'src/examples/qc/response_h2o_tpa_seeded'
work=Path(tempfile.mkdtemp(prefix='shared-seed-'))
print('artifacts:',work,flush=True)
deck=(case/'response_h2o_tpa_seeded.in').read_text()
deck=re.sub(r'convergence.criteria\s+\[[^\n]*\]','convergence.criteria [bsh_residual]',deck)
deck=re.sub(r'convergence.bsh_residual\s+\S+','convergence.bsh_residual 1000',deck)
deck=re.sub(r'(?m)^\s*dconv\s+\S+','    convergence.criteria [bsh_residual]\n    convergence.bsh_residual 1000',deck)
deck=re.sub(r'(?m)^(\s*maxiter)\s+\d+',r'\1 0',deck)
deck=deck.replace('excited.num_states    2','excited.num_states    1').replace('excited.tpa           true','excited.tpa           false')
deck=deck.replace('dalton.dir  dalton','dalton.dir  '+str(case/'dalton'))
(work/'input.in').write_text(deck)
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
cmd=[str(Path(a.binary).resolve()),'input.in','--wf=response','--prefix=run']
root=work/'run'
def run(name):
 r=subprocess.run(cmd,cwd=work,env=env,capture_output=True,text=True,timeout=600)
 (work/f'{name}.log').write_text(r.stdout+r.stderr)
 assert r.returncode==0,(r.stdout[-3000:],r.stderr)
 return r.stdout
def saved():
 path=root/'checkpoints/response/tasks/es:full'
 generation=json.loads((path/'latest').read_text())['generation']
 return json.loads((path/generation/'manifest.json').read_text())
first_log=run('initial')
assert 'imported excited bundle; fresh convergence checks' in first_log
first=saved()
assert first['task']['definition']['guess']['kind']=='dalton'
assert first['task']['iteration']['updates_this_attempt']==0
residual=next(m for m in first['task']['iteration']['metrics'] if m['metric']=='bsh_residual')
assert residual['status']=='measured' and len(residual['evidence']['components'])==10
run('reuse')
assert saved()==first
attempt=json.loads((root/'timings.json').read_text())['attempts'][-1]
def nodes(n):
 yield n
 for c in n.get('children',[]): yield from nodes(c)
assert not any(n['name'] in ['import_guesses','imported_guess','iteration'] for n in nodes(attempt['timing']))
print('shared DALTON guess/reuse checks passed',flush=True)
