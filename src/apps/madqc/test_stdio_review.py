#!/usr/bin/env python3
"""Retired response controls must reject through the real JSON input entry point."""
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
a = p.parse_args()
binary = str(Path(a.binary).resolve())
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
work = Path(tempfile.mkdtemp(prefix='stdio-json-review-'))
print('artifacts:', work, flush=True)
for key, value in {'stall.window': 3, 'stall.ratio': 0.2, 'seed.start_rung': 'fine',
                   'seed.es_warmup': True, 'excited.maxiter': 2,
                   'progress.window': 3, 'progress.minimum_relative_improvement': 0.1,
                   'accuracy.bsh_residual': 1e-4, 'accuracy.kinetic_trace_change': 1e-4,
                   'convergence.use_default_criteria': False,
                   'diagnostic.iterations_per_rung': 1}.items():
    (work/'input.json').write_text(json.dumps({'response': {key: value}}))
    result = subprocess.run(launch + [binary, 'input.json', '--wf=response', '--prefix=run'],
                            cwd=work, env=env, capture_output=True, text=True, timeout=60)
    (work/(key+'.log')).write_text(result.stdout+result.stderr)
    assert result.returncode != 0 and 'unknown input: '+key in result.stdout+result.stderr, key
    assert not (work/'run').exists(), 'invalid JSON must fail before calculations create output'
print('JSON response retirement passed', flush=True)
