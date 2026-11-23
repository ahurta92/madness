#!/usr/bin/env python3
"""Ordinary restart rejects model changes and continues paired response after SCF refinement.

Loose and impossible gates distinguish lifecycle behavior; existing qctest decks
supply the independent numerical references. Uses actual native archives.
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
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
work = Path(tempfile.mkdtemp(prefix='response-refinement-'))
print('artifacts:', work, flush=True)
env = dict(os.environ, MAD_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
binary = str(Path(a.binary).resolve())
root = work/'run'
store = root/'checkpoints'
counter = 0


def run(label, *, ground='', response='', molecule='', success=True):
    global counter
    (work/'input.in').write_text(f"""dft
 xc hf
 localize canon
 protocol [1e-4, 1e-6]
 convergence.bsh_residual 1000
 convergence.density_change 1000
 convergence.energy_change 1000
 maxiter 1
 print_level 0
 {ground}
end
response
 dipole.directions z
 dipole.frequencies [0.0, 0.05]
 convergence.bsh_residual 1000
 convergence.density_change 1000
 maxiter 1
 print_level 1
 {response}
end
molecule
 units atomic
 no_orient true
 He 0.1 0 0
 {molecule}
end
""")
    r = subprocess.run(launch+[binary, 'input.in', '--wf=response', '--prefix=run'],
                       cwd=work, env=env, capture_output=True, text=True, timeout=240)
    counter += 1
    text = r.stdout+r.stderr
    (work/f'{counter:02d}-{label}.log').write_text(text)
    assert (r.returncode == 0) == success, (label, r.returncode, text[-6000:])
    return text


def selected(calculation, task):
    path = store/calculation/'tasks'/task
    generation = json.loads((path/'latest').read_text())['generation']
    return json.loads((path/generation/'manifest.json').read_text())


def states():
    return [selected('response', f'fd:dipole_z@f{freq:.5f}') for freq in (0., .05)]


def metric(record, name):
    return next(m for m in record['task']['iteration']['metrics'] if m['metric'] == name)


def sources(records):
    return [r['task']['inputs'] for r in records]


def payloads():
    return {str(f.relative_to(store)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in store.rglob('*.00000') if f.is_file()}


run('initial')
initial = states()
assert all(r['task']['definition']['protocol'] == {'threshold': 1e-6, 'k': 8} for r in initial)
# A real SCF update fails its impossible request, leaving dependent response untouched.
run('refine-ground', ground='convergence.bsh_residual 0\n maxiter 1', success=False)
assert states() == initial
updated_ground = selected('ground_state', 'scf')
assert updated_ground['task']['iteration']['updates_this_attempt'] == 1
# The updated SCF now qualifies the relaxed request; its dependent states must be checked anew.
out = run('continue-response')
continued = states()
assert sources(continued) != sources(initial)
for previous, current in zip(initial, continued):
    assert current['task']['inputs'][0]['generation'] == updated_ground['generation']
    assert current['task']['definition']['protocol'] == previous['task']['definition']['protocol']
    # Checked anew with one update: change metrics need a predecessor to meet loose targets.
    assert current['task']['iteration']['updates_this_attempt'] == 1
    assert metric(current, 'density_change')['status'] == 'measured'
    assert metric(current, 'kinetic_change')['status'] == 'measured'
    assert current['task']['outcome']['complete']
# Existing measurements may qualify another request without a numerical update.
run('reuse', response='convergence.bsh_residual 500')
assert states() == continued
# Numerical smoothing changes are allowed, with fresh measurements and source binding.
run('smoothing', molecule='eprec 1e-8')
smoothed = states()
assert sources(smoothed) != sources(continued)
assert all(metric(r, 'density_change')['status'] == 'measured' and r['task']['iteration']['updates_this_attempt'] == 1 for r in smoothed)
# Same model, changed response request: native functions survive, comparison history resets.
run('tight-response', molecule='eprec 1e-8', response='convergence.bsh_residual 0\n maxiter 0', success=False)
tightened = states()
for previous, current in zip(smoothed, tightened):
    assert not current['task']['outcome']['complete']
    assert metric(current, 'density_change')['status'] == 'unavailable'
    assert abs(metric(current, 'bsh_residual')['value']-metric(previous, 'bsh_residual')['value']) < 1e-10
results = json.loads((root/'results.json').read_text())
response_result = next(c for c in results['calculations'] if c['calculation_id'] == 'response')
assert all(t['status'] != 'complete' for t in response_result['properties']['polarizability']['values'])
# Model rejections must not overwrite native payloads or response provenance.
for label, ground, molecule in [('localization', 'localize boys', ''),
                                ('functional', 'xc lda', ''),
                                ('cell', 'l 30', ''),
                                ('field', '', 'field [0,0,0.01]')]:
    before_payloads = payloads()
    before_states = states()
    output = run(label, ground=ground, molecule='eprec 1e-8\n'+molecule, success=False)
    assert 'model' in output or 'conversion' in output, output[-4000:]
    assert payloads() == before_payloads and states() == before_states
print('same-model SCF/response refinement and model rejection passed', flush=True)
