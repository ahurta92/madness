#!/usr/bin/env python3
"""Two-photon absorption on a named molresponse stage (D-113): H2, one RPA
root (two-photon dark by symmetry, so the test checks the plumbing, the
record's shape and its identities, not a value). requested_properties [excited_states, two_photon] adds the dipole legs
at half the root's frequency (the engine's derived legs); the record per root
holds the transition tensor S, D_f, D_g, D_linear, D_circular and R (the
theory note two-photon-normalization-sign), the root's own convergence in the
bundle and its three legs. The shared verdict judges the bundle and every leg;
--continue loads them all; a TDA bundle with two_photon is refused."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)
import named_workdir  # noqa: E402  (D-108: removed when the test passes)

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = named_workdir.work_directory('named-two-photon-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]

base = (Path(a.source) / 'named_h2_es_rpa.in').read_text()
deck = base.split('molresponse excited')[0] + """molresponse excited
  reference gs
  protocol [1e-4]
  dconv 1e-4
  maxiter 20
  requested_properties [excited_states, two_photon]
  excited_states.roots 1
  excited_states.method rpa
end
"""


def run(label, root, text, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=1200)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


root = work / 'H2'
run('two-photon', root, deck)
stage = named_records.stage_record(root / 'excited')
assert stage['status'] == 'completed', stage['status']
tp = stage['properties']['two_photon']
assert tp['converged'] is True and tp['units'] == 'au', tp
assert 'two-photon-normalization-sign' in tp['convention'], tp['convention']
assert len(tp['roots']) == 1, tp['roots']
r = tp['roots'][0]
for key in ('root', 'omega', 'photon_frequency', 'tensor', 'D_f', 'D_g', 'D_linear',
            'D_circular', 'R', 'root_convergence', 'legs', 'legs_converged',
            'max_leg_bsh_residual'):
    assert key in r, (key, sorted(r))
assert not any('sigma' in key for key in r), sorted(r)  # D-113 (b): no GM cross sections
assert r['photon_frequency'] == 0.5 * r['omega'], (r['photon_frequency'], r['omega'])
assert abs(r['D_linear'] - (2 * r['D_f'] + 4 * r['D_g'])) <= 1e-12 * max(1.0, abs(r['D_linear'])), r
assert abs(r['D_circular'] - (-2 * r['D_f'] + 6 * r['D_g'])) <= 1e-12 * max(1.0, abs(r['D_circular'])), r
S = r['tensor']
assert len(S) == 3 and all(len(row) == 3 for row in S), S
trace = S[0][0] + S[1][1] + S[2][2]
assert abs(r['D_f'] - trace * trace / 30) <= 1e-10 * max(1.0, abs(r['D_f'])), (r['D_f'], trace)
assert r['root_convergence']['converged'] is True, r['root_convergence']
assert r['root_convergence']['density_residual'] is not None
assert r['root_convergence']['excitation_energy_change'] is not None
assert len(r['legs']) == 3 and r['legs_converged'] is True, r['legs']

# The shared verdict judged the bundle and the three legs at half its omega
# (requested_accuracy is keyed by the states it judged); the stage completed,
# so every one of them converged.
judged = list(stage['convergence_report']['requested_accuracy'])
legs = [k for k in judged if k.startswith('dipole_')]
assert len(legs) == 3 and any(k.startswith('es:full_n1') for k in judged), judged
half = r['photon_frequency']
assert all(k.endswith('|f%.5f' % half) for k in legs), (legs, half)

# --continue loads the bundle and the legs: no iterations.
log = run('two-photon-continue', root, deck, resume=True)
again = named_records.stage_record(root / 'excited')
assert again['work'] == 'loaded' and again['iterations_this_run'] == 0, (again['work'], again['iterations_this_run'])
assert again['continuation']['planner'] == 'request_complete', again['continuation']
assert again['properties']['two_photon']['roots'][0]['D_linear'] == r['D_linear']

# A request that adds a static alpha is planned state by state: the loaded
# bundle keeps its omegas, so its three legs at half of it load too (the
# planner lists them from the saved bundle) and only the new leg is solved.
more = deck.replace('requested_properties [excited_states, two_photon]',
                    'requested_properties [excited_states, two_photon, polarizability]'
                    '\n  polarizability.frequencies [0.0]\n  polarizability.frequency_units au'
                    '\n  polarizability.components [xx]')
assert more != deck
run('two-photon-alpha', root, more, resume=True)
grown = named_records.stage_record(root / 'excited')
plans = {s['state']: s['action'] for s in grown['continuation']['states']}
half_legs = {s: a for s, a in plans.items() if s.startswith('dipole_') and not s.endswith('|0')}
assert len(half_legs) == 3 and all(a == 'loaded' for a in half_legs.values()), plans
# The new static leg starts from the nearest saved dipole_x frequency (R4).
assert plans.get('es:full_n1') == 'loaded' and plans.get('dipole_x|0') in ('seeded', 'fresh'), plans
assert grown['status'] == 'completed', grown['status']
# The engine assembles the value again from the loaded bundle and legs. H2's
# lowest root is two-photon dark by symmetry (a Sigma_u+ state of a
# centrosymmetric molecule), so D is numerical zero (~1e-24) and the two
# assemblies agree only to that noise; the water twins check values (D-113).
assert abs(grown['properties']['two_photon']['roots'][0]['D_linear'] - r['D_linear']) < 1e-20

# A TDA bundle with two_photon is refused at the deck check (D-113 (a)).
tda = deck.replace('excited_states.method rpa', 'excited_states.method tda')
log = run('two-photon-tda', work / 'Tda', tda, success=False)
assert 'two_photon needs excited_states in the same stage with excited_states.method rpa' in log, log[-3000:]

print(json.dumps({'omega': r['omega'], 'D_linear': r['D_linear'], 'D_circular': r['D_circular'],
                  'mpi_ranks': a.ranks}, indent=2))
print('PASS')
