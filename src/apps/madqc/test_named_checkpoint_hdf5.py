#!/usr/bin/env python3
"""A response HDF5 mirror must not add to or overwrite a published generation."""
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
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='checkpoint-hdf5-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MADRESPONSE_IO_HDF5='1')
text = (Path(a.source) / 'named_h2.in').read_text()
text = text.replace('[1e-4, 1e-6]', '[1e-4]').replace('dconv 1e-5', 'dconv 1e-3')
text = text.replace('dconv 1e-4', 'dconv 1e-3').replace('[0.0, 0.02]', '[0.0]')
root = work / 'Run'
def run(label, deck, resume=False):
    inp = work / (label + '.in')
    inp.write_text(deck)
    with (work / (label + '.log')).open('w') as log:
        proc = subprocess.run([binary, str(inp), '--prefix=' + str(root)] +
                              (['--continue'] if resume else []), cwd=work, env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=300)
    assert proc.returncode == 0, (work / (label + '.log')).read_text()[-4000:]
run('ground', text.split('molresponse optical')[0])
store = root / 'gs/checkpoints'
generation = json.loads((store / 'latest').read_text())['generation']
directory = store / generation
def hashes():
    return {str(f.relative_to(directory)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in directory.rglob('*') if f.is_file()}
before = hashes()
response_generation = None
for label in ('response', 'repeat', 'additional-frequency'):
    requested = text.replace('[0.0]', '[0.0, 0.02]') if label == 'additional-frequency' else text
    run(label, requested, True)
    assert hashes() == before, 'response HDF5 mirror changed an immutable generation'
    response = json.loads((root / 'optical/results.json').read_text())
    assert response['reference']['state_id'] == generation
    if label == 'repeat':
        assert response['work'] == 'loaded' and response['iterations_this_run'] == 0
        assert response['checkpoint_generation'] == response_generation
    else:
        assert list(Path(response['output_directory']).rglob('*.h5')), 'HDF5 path was not exercised'
        assert response['iterations_this_run'] > 0
    if label == 'response':
        response_generation = response['checkpoint_generation']
    response_store = root / 'optical/checkpoints'
    selected = response_store / response['checkpoint_generation']
    manifest = json.loads((selected / 'manifest.json').read_text())
    assert all(e['backend'] == 'native' for c in manifest['numerical']['engine_metadata']['fd_states'].values()
               for rung in c.values() for e in rung.values())
print('HDF5 response/repeat preserved the selected generation', flush=True)
