#!/usr/bin/env python3
"""Multiple independent states, mixed progress, partial output and a real interrupted ladder."""
from response_test_checks import assert_same_alpha_properties
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time

p=argparse.ArgumentParser()
p.add_argument('--binary',required=True);p.add_argument('--source',required=True)
p.add_argument('--wiring-only',action='store_true');p.add_argument('--ranks',type=int,default=1);p.add_argument('--mpiexec',default='mpiexec')
a=p.parse_args();work=Path(tempfile.mkdtemp(prefix='response-ladder-'));print('artifacts:',work,flush=True)
source=Path(a.source).resolve()
(work/'input.in').write_text((source/'src/examples/qc/response_shared_ladder/response_shared_ladder.in').read_text())
root=work/'run';tasks=root/'checkpoints/response/tasks'
env=dict(os.environ,MAD_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
launch=[] if a.ranks==1 else [a.mpiexec,'-n',str(a.ranks)]
command=launch+[str(Path(a.binary).resolve()),'input.in','--wf=response','--prefix=run']
count=0

def task(frequency,axis='z'):return f'fd:dipole_{axis}@f{frequency:.5f}'
def saved(identity):
    g=json.loads((tasks/identity/'latest').read_text())['generation']
    return g,json.loads((tasks/identity/g/'manifest.json').read_text())
def residual(record):return next(m['value'] for m in record['task']['iteration']['metrics'] if m['metric']=='bsh_residual')
def response_result():return next(c for c in json.loads((root/'results.json').read_text())['calculations'] if c['calculation_id']=='response')
def history():
    return [json.loads(line) for f in sorted((work/'run/analysis/response').glob('*.jsonl'),key=lambda f:f.stat().st_mtime_ns) for line in f.read_text().splitlines()]
def run(options='',success=True,interrupt=None):
    global count
    if a.wiring_only and 'convergence.bsh_residual=' not in options:
        options=loose+(';' + options if options else '')
    count+=1
    with (work/f'{count:02d}.log').open('w') as log:
        proc=subprocess.Popen(command+(['--response='+options] if options else []),cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        protected=None
        try:
            if interrupt:
                deadline=time.monotonic()+240
                while proc.poll() is None and time.monotonic()<deadline:
                    pointer=tasks/interrupt/'latest'
                    if pointer.exists() and json.loads(pointer.read_text())['generation']:
                        _,m=saved(interrupt)
                        if m['task']['definition']['protocol']['threshold']==1e-5:
                            protected=tasks/interrupt;protected.chmod(0o500);break
                    time.sleep(0.01)
                assert protected,'did not observe intermediate checkpoint'
            code=proc.wait(timeout=600)
        finally:
            if protected:protected.chmod(0o700)
            if proc.poll() is None:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
    assert (code==0)==success,(code,(work/f'{count:02d}.log').read_text()[-6000:])

# The deck also writes a kinetic target, so it is a criterion and gets the loose target too.
loose='convergence.bsh_residual=1000;convergence.density_change=1000;convergence.kinetic_change=1000;maxiter=1'
run(loose+';dipole.frequencies=[0.05]')
dynamic=task(.05);static=task(0)
original,_=saved(dynamic)
# One update per level: density change has a predecessor only after it.
assert all((next(m for m in r['metrics'] if m['metric']=='density_change')['status']=='unavailable')==(r['iteration']==0) for r in history() if r['task_id']==dynamic)
assert saved(dynamic)[1]['task']['definition']['protocol']=={'threshold':1e-6,'k':8}
assert [r['protocol']['threshold'] for r in history() if r['task_id']==dynamic]==[1e-4,1e-4,1e-5,1e-5,1e-6,1e-6]
run(loose)
assert saved(dynamic)[0]==original,'new state restarted an already finished state'
assert [r['protocol']['threshold'] for r in history() if r['task_id']==static]==[1e-4,1e-4,1e-5,1e-5,1e-6,1e-6]
assert len(list(tasks.glob('*/*/*.00000')))==2,'only latest payload per state should remain'
before={t:saved(t)[0] for t in [static,dynamic]}
values={t:residual(saved(t)[1]) for t in before};gate=sum(values.values())/2
assert min(values.values())<gate<max(values.values())
run(f'convergence.bsh_residual={gate:.17g};convergence.density_change=1000;convergence.kinetic_change=1000;maxiter=0',success=False)
for t,value in values.items():assert (saved(t)[0]==before[t])==(value<=gate)
rows=response_result()['properties']['polarizability']['values']
assert sum('zz' in r['components'] for r in rows)==1
assert sum('zz' in r['availability'] for r in rows)==1
run()
converged={t:saved(t)[0] for t in before}
if not a.wiring_only:
    reference=json.loads((source/'src/examples/qc/response_he_alpha/reference/response_he_alpha.calc_info.json').read_text())['tasks'][1]['properties']['response_properties']['alpha']['1e-04_k6']
    for row in response_result()['properties']['polarizability']['values']:
        expected=next(r['alpha'][0][0] for r in reference if r['omega']==row['frequency'])
        assert abs(row['components']['zz']-expected)<1e-3
        print('alpha_zz',row['frequency'],row['components']['zz'],'reference',expected,flush=True)
digests={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in tasks.glob('*/*/*.00000')}
rows_before=len(history())
before_reuse=response_result()
original_sources={t:(saved(t)[0],saved(t)[1]['attempt_id']) for t in before}
# The native manifests are authoritative even when analysis metadata is lost.
(work/'run/analysis/response/response_metadata.json').unlink()
run()
after_reuse=response_result()
assert after_reuse['attempt_id']!=before_reuse['attempt_id']
assert_same_alpha_properties(before_reuse['properties'],after_reuse['properties'])
sources=json.loads((root/'checkpoints/response/attempts'/after_reuse['attempt_id']/'result_sources.json').read_text())['states']
assert set(sources)==set(original_sources)
for t,(generation,attempt) in original_sources.items():
    assert (sources[t]['generation'],sources[t]['attempt_id'])==(generation,attempt)
assert {t:saved(t)[0] for t in before}==converged and len(history())==rows_before
assert all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in digests.items())
added=task(.1)
run('dipole.frequencies=[0,0.05,0.1]',success=False,interrupt=added)
assert saved(added)[1]['task']['iteration']['protocol_complete']
assert saved(added)[1]['task']['definition']['protocol']=={'threshold':1e-5,'k':8}
rows_before=len(history());run('dipole.frequencies=[0,0.05,0.1]')
new_rows=history()[rows_before:]
assert new_rows and all(r['task_id']==added and r['protocol']['threshold']==1e-6 for r in new_rows)
assert {t:saved(t)[0] for t in before}==converged
run('dipole.directions=xyz')
result=response_result();rows=result['properties']['polarizability']['values']
assert len(rows)==2 and all(len(r['components'])==9 and not r['availability'] for r in rows)
assert sum(len(group['states']) for group in result['state_metrics'])==6
assert {t:saved(t)[0] for t in before}==converged
if not a.wiring_only:
    for row in rows:
        for axis in 'xy':assert abs(row['components'][axis+axis]-row['components']['zz'])<1e-3
        for component,value in row['components'].items():
            if component[0]!=component[1]:assert abs(value)<1e-3
print('multi-state protocol/restart and grouped tensor checks passed')
