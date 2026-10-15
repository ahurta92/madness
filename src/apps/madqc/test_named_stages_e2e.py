#!/usr/bin/env python3
"""Fresh named-stage/legacy numerical parity and output/identity isolation."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True)
parser.add_argument('--source', required=True)
parser.add_argument('--ranks', type=int, default=1)
parser.add_argument('--mpiexec', default='mpiexec')
args = parser.parse_args()
binary = str(Path(args.binary).resolve())
source = Path(args.source)
work = Path(tempfile.mkdtemp(prefix='named-stages-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ)
env.setdefault('MAD_NUM_THREADS', '2')
env.setdefault('OPENBLAS_NUM_THREADS', '1')
env.setdefault('OMP_NUM_THREADS', '1')
launch = [] if args.ranks == 1 else [args.mpiexec, '-n', str(args.ranks)]

def run(name, deck, options, success=True):
    directory = work / name
    directory.mkdir()
    (directory / 'deck.in').write_text(deck)
    start = time.monotonic()
    with (directory / 'run.log').open('w') as log:
        proc = subprocess.run(launch + [binary, 'deck.in', *options], cwd=directory,
                              env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600)
    print(name, 'rc', proc.returncode, 'seconds', round(time.monotonic()-start, 2), flush=True)
    assert (proc.returncode == 0) == success, (directory / 'run.log').read_text()[-6000:]
    return directory

named_deck = (source / 'named_h2.in').read_text()
# Two response instances exercise distinct requests/targets on one source.
second = '''\nmolresponse static_check
 reference gs
 protocol [1e-4,1e-6]
 dconv 2e-4
 maxiter 30
 kain true
 requested_properties [polarizability]
 polarizability.frequencies [0.0]
 polarizability.frequency_units au
 polarizability.components [zz]
end
'''
named = run('named', named_deck + second, ['--prefix=CaseRun'])
assert sorted(p.name for p in named.iterdir()) == ['CaseRun', 'deck.in', 'run.log'], 'stage output escaped its root'
results = json.loads((named / 'CaseRun/results.json').read_text())
gs, response, second_response = results['stages']
assert [s['stage'] for s in results['stages']] == ['gs', 'optical', 'static_check']
assert all(s['status'] == 'completed' for s in results['stages'])
assert len({s['state_id'] for s in results['stages']}) == 3
for s in results['stages']:
    assert (named / 'CaseRun' / s['stage'] / 'results.json').exists()
assert gs['method'] == 'HF' and response['method'] == 'TDHF'
assert gs['properties']['energy']['units'] == 'hartree'
assert gs['properties']['dipole']['units'] == 'e*a0'
assert max(abs(v) for v in gs['properties']['dipole']['value']) < 1e-8, 'H2 dipole must vanish by symmetry'
assert response['reference']['state_id'] == gs['state_id']
assert second_response['reference']['state_id'] == gs['state_id']
assert response['targets']['dconv'] == 1e-4 and second_response['targets']['dconv'] == 2e-4
alpha = response['properties']['polarizability']['values']
assert [v['frequency'] for v in alpha] == [0.0, 0.02]
assert all(v['units'] == 'a0^3' and set(v['components']) == {'zz'} for v in alpha)
assert len(second_response['properties']['polarizability']['values']) == 1

legacy_deck = '''dft
 xc hf
 L 20
 localize canon
 protocol [1e-4,1e-6]
 dconv 1e-5
 econv 1e-6
 maxiter 30
 dipole true
end
response
 dconv 1e-4
 maxiter 30
 kain true
 requested_properties [polarizability]
 dipole.frequencies [0.0,0.02]
 dipole.directions z
end
molecule
 units atomic
 eprec 1e-6
 no_orient true
 H 0 0 -0.7
 H 0 0 0.7
end
'''
legacy = run('legacy', legacy_deck, ['--wf=response', '--prefix=comparison'])
old = json.loads((legacy / 'comparison.calc_info.json').read_text())['tasks']
assert abs(gs['properties']['energy']['value'] - old[0]['properties']['energy']) < 1e-9
old_alpha = old[1]['properties']['response_properties']['alpha']['1e-06_k8']
for a in alpha:
    b = next(r for r in old_alpha if abs(r['omega'] - a['frequency']) < 1e-12)
    assert abs(a['components']['zz'] - b['alpha'][0][0]) < 1e-6
# Existing output roots cannot be accidentally reused during the fresh-only slice.
with (named / 'repeat.log').open('w') as log:
    repeat = subprocess.run(launch + [binary, 'deck.in', '--prefix=CaseRun'], cwd=named,
                            env=env, stdout=log, stderr=subprocess.STDOUT, timeout=60)
assert repeat.returncode != 0
assert 'new output directory' in (named / 'repeat.log').read_text()
assert json.loads((named / 'CaseRun/results.json').read_text()) == results
bad = run('duplicate', named_deck + second.replace('static_check', 'gs'), ['--prefix=bad'], False)
assert 'duplicate name' in (bad / 'run.log').read_text()
assert not (bad / 'bad').exists(), 'invalid deck began a calculation'
# An omitted L is derived again after SCF centers the geometry. Response must
# bind that actual cell, not the larger pre-centering estimate from the parser.
shifted_deck = named_deck.replace('  no_orient true\n', '').replace('  L 20\n', '')
shifted_deck = shifted_deck.replace('H 0 0 -0.7', 'H 0 0 9.3').replace('H 0 0 0.7', 'H 0 0 10.7')
shifted_deck = shifted_deck.replace('[1e-4, 1e-6]', '[1e-4]').replace('dconv 1e-5', 'dconv 1e-3')
shifted_deck = shifted_deck.replace('dconv 1e-4', 'dconv 1e-3').replace('[0.0, 0.02]', '[0.0]')
shifted = run('shifted', shifted_deck, ['--prefix=Shifted'])
shifted_results = json.loads((shifted / 'Shifted/results.json').read_text())['stages']
assert abs(shifted_results[0]['representation']['L'] - 50.7) < 1e-8
assert shifted_results[1]['reference_representation']['L'] == shifted_results[0]['representation']['L']
protocol_lines = [line for line in (shifted / 'run.log').read_text().splitlines() if line.startswith('PROTOCOL_SET')]
assert protocol_lines and all(abs(float(line.split('L=')[1].split()[0])-50.7) < 1e-8 for line in protocol_lines)
print(json.dumps({'energy': gs['properties']['energy'], 'dipole': gs['properties']['dipole'],
                  'alpha': alpha, 'mpi_ranks': args.ranks}, indent=2))
