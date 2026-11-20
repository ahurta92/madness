#!/usr/bin/env python3
"""SCF and response retain measured current-invocation timings in one file."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--ranks', type=int, default=1)
parser.add_argument('--mpiexec', default='mpiexec')
args = parser.parse_args()
work = Path(tempfile.mkdtemp(prefix='shared-workflow-timings-'))
print('artifacts:', work, flush=True)
(work / 'input.in').write_text('''dft
 xc hf
 localize canon
 protocol [1e-4]
 l 30
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 1
end
response
 dipole true
 dipole.directions z
 dipole.frequencies [0.0]
 maxiter 0
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 print_level 1
end
molecule
 units atomic
 He 0 0 0
end
''')
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]
root = work / 'run'
def run(name, extra=(), success=True):
    result = subprocess.run(launch + [str(Path(args.binary).resolve()), 'input.in', '--wf=response',
                            '--prefix=run', *extra], cwd=work, env=env, capture_output=True,
                            text=True, timeout=180)
    text = result.stdout + result.stderr
    (work / (name + '.log')).write_text(text)
    assert (result.returncode == 0) == success, text[-3000:]
    return text

def walk(node):
    yield node
    for child in node.get('children', []):
        yield from walk(child)

def records():
    return json.loads((root / 'timings.json').read_text())['attempts']

def calculation(record):
    return record['timing']['children'][0]['calculation_id']

def iterations(record):
    return [n for n in walk(record['timing']) if n['name'] == 'iteration']

def validate(record):
    assert record['execution']['mpi_ranks'] == args.ranks
    assert record['timing']['clock_scope'] == 'driver_rank'
    for n in walk(record['timing']):
        assert math.isfinite(n['seconds']) and n['seconds'] >= 0, n
        for child in n.get('children', []):
            assert child['seconds'] <= n['seconds'] + 1e-8, (n, child)

def compare_console(text, record):
    rows = re.findall(r'Timing SCF:.*?\(iteration (\d+)\)\nIteration wall time \(seconds, driver rank\)\n(.*?)(?=\n[^ ]|\Z)', text, re.S)
    saved = iterations(record)
    assert len(rows) == len(saved), (rows, saved)
    for (index, body), node in zip(rows, saved):
        assert int(index) == node['iteration']
        printed = re.findall(r'^  (.+?)\s+(\d+\.\d+) *$', body, re.M)
        assert len(printed) == len(node['children']) + 1, (printed, node)
        for (_, value), measured in zip(printed, [c['seconds'] for c in node['children']] + [node['seconds']]):
            assert abs(float(value) - measured) <= 0.000501, (value, measured)

first_output = run('fresh')
first = records()
assert [calculation(r) for r in first] == ['ground_state', 'response'], 'shared timings must include SCF before response'
for record in first:
    validate(record)
    assert len(iterations(record)) == 1
    assert 'advance' not in [c['name'] for c in iterations(record)[0]['children']]
compare_console(first_output, first[0])
sources = {str(p.relative_to(root)): p.read_bytes() for p in (root / 'checkpoints').rglob('latest')}
run('reuse')
reused = records()
assert len(reused) == 4 and reused[:2] == first
for record in reused[2:]:
    validate(record)
    assert not iterations(record), 'native reuse must not report iterations'
assert sources == {str(p.relative_to(root)): p.read_bytes() for p in (root / 'checkpoints').rglob('latest')}
last_output = run('tighten', ['--dft=convergence.bsh_residual=0;maxiter=1'], success=False)
last = records()
assert len(last) == 5 and last[:4] == reused
validate(last[-1])
measured = iterations(last[-1])
assert len(measured) == 2 and [n['iteration'] for n in measured] == [0, 1]
assert 'advance' in [c['name'] for c in measured[0]['children']]
assert 'advance' not in [c['name'] for c in measured[-1]['children']]
compare_console(last_output, last[-1])
print('shared workflow timing checks passed', flush=True)
