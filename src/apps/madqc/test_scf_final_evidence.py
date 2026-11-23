#!/usr/bin/env python3
"""Endpoint/result metrics identify exactly the native generation SCF saved."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--binary',required=True)
p.add_argument('--source')  # kept for the existing test registration
args=p.parse_args()
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS','2')
with tempfile.TemporaryDirectory(prefix='scf-final-evidence-') as scratch:
    work=Path(scratch)
    (work/'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 maxiter 2
 convergence.bsh_residual 1000
 convergence.energy_change 0
 convergence.density_change 1000
 convergence.kinetic_change 1000
 print_level 0
end
molecule
 units atomic
 He 0 0 0
end
''')
    run=subprocess.run([str(Path(args.binary).resolve()),'input.in','--prefix=run'],cwd=work,env=env,capture_output=True,text=True,timeout=90)
    assert run.returncode==1,run.stdout+run.stderr
    root=work/'run'
    task=root/'checkpoints/ground_state/tasks/scf'
    generation=json.loads((task/'latest').read_text())['generation']
    manifest=json.loads((task/generation/'manifest.json').read_text())
    saved=manifest['task']['iteration']
    assert saved['updates_this_attempt']==2 and not saved['protocol_complete']
    result=json.loads((root/'results.json').read_text())['calculations'][0]
    energy=result['properties']['total_energy']
    assert energy['status']=='unavailable' and 'value' not in energy
    reported=result['state_metrics'][0]['states']['ground']
    for metric in saved['metrics']:
        assert metric['status']=='measured',metric
        assert reported['values'][metric['metric']]==metric['value'],metric
    assert (task/generation/'scf.restartdata.00000').exists()
print('SCF final evidence preserved')
