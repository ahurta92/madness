#!/usr/bin/env python3
"""D-061: econv on a molresponse stage that requests excited_states is the
target of the excitation energy change per root (D8: a written value only),
reported under the engine check and never enforced (B8): the solve is the
same with and without it, and on a stage without excited_states it is refused.
"""
import argparse
import json
import sys
import os
from pathlib import Path
import subprocess
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
work = named_workdir.work_directory('named-es-econv-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
without_econv = (Path(a.source) / 'named_h2_es_tda.in').read_text()
assert 'econv' not in without_econv.split('molresponse excited')[1]
with_econv = without_econv.replace('  excited_states.method tda\n',
                                   '  excited_states.method tda\n  econv 1e-6\n')
assert with_econv != without_econv


def run(label, root, text, expect=0):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)],
                           cwd=work, env=env, stdout=out, stderr=subprocess.STDOUT,
                           timeout=900)
    if expect is not None:
        assert (r.returncode == 0) == (expect == 0), log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


def excited(root):
    stages = {s['stage']: s for s in named_records.stages(root)}
    return stages['excited']


# The written econv is reported, not enforced: the engine check still stops
# the bundle, so the solve is the same as without it (B4).
fresh = work / 'WithEconv'
run('with-econv', fresh, with_econv)
plain = work / 'WithoutEconv'
run('without-econv', plain, without_econv)
es_with, es_without = excited(fresh)['properties']['excited_states'], excited(plain)['properties']['excited_states']
assert es_with['converged'] and es_without['converged'], (es_with, es_without)
assert es_with['iterations_at_protocol_step'] == es_without['iterations_at_protocol_step'], (es_with, es_without)
assert excited(fresh)['iterations_this_run'] == excited(plain)['iterations_this_run']
for root_with, root_without in zip(es_with['roots'], es_without['roots']):
    assert root_with['root_id'] == root_without['root_id']
    assert abs(root_with['omega'] - root_without['omega']) < 1e-10, (root_with, root_without)
    # The last step's eigenvalue change is a difference of two omegas, so it
    # carries their roundoff; the same tolerance as omega itself.
    assert abs(root_with['excitation_energy_change'] - root_without['excitation_energy_change']) < 1e-10, (root_with, root_without)
    assert root_with['excitation_energy_change'] <= 1e-6, root_with  # H2 TDA ends far inside econv 1e-6

# D-096: the engine's screen at the default level (print_level unwritten,
# Normal): the first iteration's density change is unavailable, not zero;
# one identified start and end line per state and rung; the save, memory
# and dispatch lines are detail at verbose.
screen = (work / 'with-econv.log').read_text()
assert 'dρ = unavailable (first sample)' in screen, 'the first sample must be marked unavailable'
assert screen.count('dipole_z f0.00000 at 1e-04_k6: start, seed') == 1, screen.count('dipole_z f0.00000 at 1e-04_k6: start')
assert screen.count('dipole_z f0.00000 at 1e-04_k6: converged in') == 1, screen.count('dipole_z f0.00000 at 1e-04_k6: converged in')
assert screen.count('es:tda_n1 at 1e-04_k6: converged in') == 1, screen.count('es:tda_n1 at 1e-04_k6: converged in')
for noise in ('[SAVE]', 'MEMORY_HWM', '[CALC] run_protocol', '[CALC] fd solve', '[KAIN]', 'Converged in'):
    assert noise not in screen, noise + ' must be verbose-only'

# The record describes the request beside the engine's gate and points at
# the shared verdict; it carries no status of its own (D-061).
block = es_with['requested_omega_target']
assert block['target'] == 1e-6 and block['set_by'] == 'user', block
assert block['units'] == 'hartree' and block['relation'] == 'less_than_or_equal', block
assert block['worst']['value'] == max(r['excitation_energy_change'] for r in es_with['roots']), block
assert block['worst']['root_id'] in {r['root_id'] for r in es_with['roots']}, block
assert all(r['excitation_energy_change_within_requested'] is True for r in es_with['roots']), es_with['roots']
assert not {'status', 'enforced', 'outcome'} & set(block), block
pointer = block['verdict_in']
assert pointer['path'] == 'convergence_report.requested_accuracy', pointer
assert pointer['criterion'] == 'excitation_energy_change', pointer
# S4d2: one requested_accuracy entry per state key on a response stage; the
# pointer names the bundle's entry (its status is the shared verdict's).
accuracy = excited(fresh)['convergence_report']['requested_accuracy']
assert pointer['state_key'] in accuracy, (pointer, sorted(accuracy))
plain_block = es_without['requested_omega_target']
assert plain_block['target'] is None and plain_block['set_by'] == 'not_set', plain_block
assert all('excitation_energy_change_within_requested' not in r for r in es_without['roots'])

# D-092: the same deck in criteria mode. The bundle's solve stops on the
# engine check in both modes, so it ends at the same iteration with the same
# omega; the written econv is the stage verdict's, judged on the bundle's
# recorded eigenvalue change (D-082): the stage converges and the pointer's
# entry under requested_accuracy is satisfied.
criteria_deck = with_econv.replace(' econv 1e-6\n', ' econv 1e-6\n convergence.use_default_criteria false\n')
assert criteria_deck != with_econv
crit = work / 'CriteriaEconv'
run('criteria-econv', crit, criteria_deck)
stage_crit = excited(crit)
assert stage_crit['convergence_report']['decided_by'] == 'requested_criteria', stage_crit['convergence_report']['decided_by']
assert stage_crit['status'] == 'completed', (stage_crit['status'], stage_crit['convergence_report'])
pointer_crit = stage_crit['properties']['excited_states']['requested_omega_target']['verdict_in']
entry_crit = stage_crit['convergence_report']['requested_accuracy'][pointer_crit['state_key']]
assert entry_crit['status'] == 'satisfied', entry_crit
es_crit = stage_crit['properties']['excited_states']
assert es_crit['iterations_at_protocol_step'] == es_with['iterations_at_protocol_step'], (es_crit['iterations_at_protocol_step'], es_with['iterations_at_protocol_step'])
assert es_crit['stop_reason'] == 'converged' and es_crit['decided_by'] == 'engine_check', (es_crit['stop_reason'], es_crit['decided_by'])
for root_crit, root_with in zip(es_crit['roots'], es_with['roots']):
    assert root_crit['root_id'] == root_with['root_id']
    assert abs(root_crit['omega'] - root_with['omega']) < 1e-10, (root_crit, root_with)
assert es_crit['requested_omega_target']['target'] == 1e-6, es_crit['requested_omega_target']

# A driven response has no local energy: econv without excited_states is
# refused at the deck check, before any stage runs.
optical = without_econv.replace('[excited_states, polarizability]', '[polarizability]')
optical = optical.replace('  excited_states.roots 1\n  excited_states.method tda\n', '  econv 1e-6\n')
assert 'excited_states' not in optical.split('molresponse excited')[1]
log = run('optical-econv', work / 'OpticalEconv', optical, expect=1)
assert 'econv is not applicable: a driven response has no local energy' in log, log[-2000:]
assert not (work / 'OpticalEconv' / 'results.json').exists()
print('PASS')
