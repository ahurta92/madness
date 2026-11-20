#!/usr/bin/env python3
"""DALTON ground seeds use the ordinary SCF driver; restart wins over reseeding."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--source',required=True)
p.add_argument('--numerical',action='store_true')
p.add_argument('--ranks',type=int,default=1)
p.add_argument('--mpiexec',default='mpiexec')
a=p.parse_args()
seed=(Path(a.source)/'src/examples/qc/response_h2o_tpa_seeded/dalton').resolve()
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS','2')
launch=[] if a.ranks==1 else [a.mpiexec,'-n',str(a.ranks)]
with tempfile.TemporaryDirectory(prefix='madqc-seed-') as scratch:
    work=Path(scratch)
    deck=f'''dft
 xc hf
 localize canon
 protocol [1e-4]
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 0
end
io
 dalton.dir {seed}
end
molecule
 eprec 1e-6
 units atomic
 no_orient true
 O 0.0000000000 0.0000005943 0.1925923500
 H 0.0000000000 1.4214240411 -0.8721230117
 H 0.0000000000 -1.4214246355 -0.8721234678
end
'''
    def run(prefix,text=deck):
        (work/'input.in').write_text(text)
        return subprocess.run(launch+[str(Path(a.binary).resolve()),'input.in','--wf=scf','--prefix='+prefix],cwd=work,env=env,capture_output=True,text=True,timeout=600)
    if a.numerical:
        accurate = deck.replace('maxiter 0','maxiter 30').replace(
            'convergence.criteria [bsh_residual]',
            'convergence.criteria [bsh_residual,density_change,energy_change]').replace(
            'convergence.bsh_residual 1000',
            'convergence.bsh_residual 5e-4\n convergence.density_change 1e-4\n convergence.energy_change 1e-6')
        cold = accurate.replace(f'io\n dalton.dir {seed}\nend\n','')
        for prefix, text in [('cold',cold),('seeded',accurate)]:
            result=run(prefix,text)
            assert result.returncode==0,result.stdout+result.stderr
        records=[]
        for prefix in ['cold','seeded']:
            root=work/prefix
            energy=json.loads((root/'results.json').read_text())['calculations'][0]['properties']['total_energy']['value']
            task=root/'checkpoints/ground_state/tasks/scf'
            generation=json.loads((task/'latest').read_text())['generation']
            iterations=json.loads((task/generation/'manifest.json').read_text())['task']['iteration']['updates_this_attempt']
            records.append((energy,iterations))
        print('cold/seeded (energy, updates):',records,flush=True)
        assert abs(records[0][0]-records[1][0])<5e-5,records
        assert records[1][1]<records[0][1],records
        print('DALTON and atomic guesses reach the same SCF energy')
        raise SystemExit(0)
    first=run('seeded')
    assert first.returncode==0,first.stdout+first.stderr
    assert '[DALTON-SEED] GS: seeding' in first.stdout, 'SCF workflow ignored DALTON seed'
    root=work/'seeded'
    pointer=root/'checkpoints/ground_state/tasks/scf/latest'
    previous=pointer.read_bytes()
    assert json.loads((root/'analysis/ground_state/seeded.gs_seed.json').read_text())['n_occ']==5
    # A native shared checkpoint is authoritative even after the import archive is removed.
    for archive in (root/'analysis/ground_state').glob('seeded.restartdata*'):
        archive.unlink()
    repeated=run('seeded')
    assert repeated.returncode==0,repeated.stdout+repeated.stderr
    assert pointer.read_bytes()==previous
    assert '[DALTON-SEED] GS: seeding' not in repeated.stdout
    moved=run('moved',deck.replace('0.1925923500','0.3925923500'))
    assert moved.returncode==1 and 'GEOMETRY FINGERPRINT MISMATCH' in moved.stdout
    assert not (work/'moved/checkpoints/ground_state/tasks/scf/latest').exists()
print('madqc seed contract passed')
