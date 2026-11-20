#!/usr/bin/env python3
"""2PA tensor identities and response reuse through the sole public driver.

This checks the engine's current records. Shared response checkpoint adoption
is a separate acceptance gate; no named-stage parser or record reader is used.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
base = (Path(a.source)/'src/examples/qc/response_h2_es_rpa/response_h2_es_rpa.in').read_text()
deck = base.replace('excited.tda            false', 'excited.tda            false\n    excited.tpa            true')
assert deck != base
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '7')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
work = Path(tempfile.mkdtemp(prefix='madqc-response-'))
print('artifacts:', work, flush=True)
(work/'input.in').write_text(deck)


def run(label):
    with (work/(label+'.log')).open('w') as log:
        proc = subprocess.run(launch+[str(Path(a.binary).resolve()),'--wf=response',
                              'input.in','--prefix=run'], cwd=work, env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=900)
    assert proc.returncode == 0, (work/(label+'.log')).read_text()[-6000:]
    return json.loads((work/'run.calc_info.json').read_text())['tasks'][1]


def tpa(task):
    rows = task['properties']['response_properties']['tpa']['1e-04_k6']
    assert len(rows) == 1
    row = rows[0]
    assert row['converged']
    assert abs(row['D_linear']-(2*row['Df']+4*row['Dg'])) < 1e-12*max(1,abs(row['D_linear']))
    assert abs(row['D_circular']-(-2*row['Df']+6*row['Dg'])) < 1e-12*max(1,abs(row['D_circular']))
    s = row['S']
    assert len(s) == 3 and all(len(x) == 3 for x in s)
    assert abs(row['Df']-sum(s[i][i] for i in range(3))**2/30) < 1e-12
    assert abs(row['Dg']-sum(x*x for r in s for x in r)/30) < 1e-12
    legs = row['row_accuracy']
    assert len(legs) == 3 and all(x['converged'] for x in legs)
    assert {x['input'] for x in legs} == {'dipole_x','dipole_y','dipole_z'}
    assert all(abs(x['freq']-0.5*row['omega']) < 1e-12 for x in legs)
    assert task['metadata']['excited_states']['1e-04_k6']['converged']
    return row


fresh = run('fresh')
first = tpa(fresh)
scf = work/'run/checkpoints/ground_state/tasks/scf/latest'
original = scf.read_bytes()
archives = {str(f.relative_to(work)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in (work/'run').rglob('*.00000')}
assert archives, 'calculation saved no native functions'
again = run('restart')
second = tpa(again)
assert scf.read_bytes() == original, 'qualified ground-state checkpoint was replaced'
assert abs(first['omega']-second['omega']) < 1e-12
assert abs(first['D_linear']-second['D_linear']) < 1e-20  # H2 root is two-photon dark
for path, digest in archives.items():
    assert hashlib.sha256((work/path).read_bytes()).hexdigest() == digest, path
assert again['metadata']['v3_diagnostics']['stop_reason'] == 'complete'
print('workflow 2PA and native restart passed')
# Preserve failed artifacts; remove successful scratch work (archives are large).
import shutil
shutil.rmtree(work)
