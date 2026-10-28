#!/usr/bin/env python3
"""Checkpoint cadence and retention of a named moldft stage (DESIGN A3), and
--continue from a kept mid-rung generation (B6).

- By default a stage publishes only at decision points (every rung end and
  stop) and keeps only those: older mid-rung generations are deleted after
  each publish. The stop manifest carries the full per-iteration history.
- checkpoint.every_iterations K also publishes every K-th iteration at a rung;
  checkpoint.keep_all true keeps every generation.
- --continue from a kept mid-rung generation resumes that rung at the next
  iteration; retention then deletes the generation it resumed from, and the
  earlier run's history moves into the new manifests.
- A deletion that fails is a warning on every rank; the stage still finishes
  (DESIGN P1: the failure is on rank 0's file work only)."""
import argparse
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='named-cadence-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
# The ground-state stage of named_h2.in: two rungs, so rung ends are tested.
deck = (Path(a.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
assert 'protocol [1e-4, 1e-6]' in deck and '  maxiter 30\n' in deck


def with_keys(*lines):
    return deck.replace('  maxiter 30\n', '  maxiter 30\n' + ''.join('  %s\n' % l for l in lines))


def run(label, root, text=deck, resume=False):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=600)
    assert r.returncode == 0, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


def result(root):
    return json.loads((root / 'gs/results.json').read_text())


def generations(root):
    """Published gs generations, oldest first; names starting with '.' (the
    publish lock, pointer temporaries, .deleting-*) are not generations."""
    store = root / 'gs/checkpoints'
    found = [json.loads((d / 'manifest.json').read_text()) for d in store.iterdir()
             if d.is_dir() and not d.name.startswith('.') and (d / 'manifest.json').exists()]
    # g-<stage>-<run start time>-<pid>-<sequence> (D-058; g-<run start
    # time>-<pid>-<sequence> before it): order by run, then by sequence. The
    # last three fields are fixed; a stage name may contain '-'.
    return sorted(found, key=lambda m: (int(m['generation'].split('-')[-3]),
                                        int(m['generation'].split('-')[-1])))


def rung(m):
    return m['position']['rung']


def is_decision_point(m):
    return (m['position']['phase'] == 'final_properties_on_saved_orbitals'
            or rung(m)['stop'] != 'in_progress')


def latest(root):
    store = root / 'gs/checkpoints'
    return json.loads((store / json.loads((store / 'latest').read_text())['generation'] /
                       'manifest.json').read_text())


# 1. Default: decision points only; nothing mid-rung kept; full history at the stop.
fresh = work / 'Default'
log = run('default', fresh)
kept = generations(fresh)
assert all(is_decision_point(m) for m in kept), [m['position'] for m in kept]
ends = [m for m in kept if m['position']['phase'] == 'after_final_rotation']
assert [rung(m)['index'] for m in ends] == [1, 2], [rung(m) for m in ends]
assert [rung(m)['stop'] for m in ends] == ['converged', 'converged'], [rung(m) for m in ends]
assert kept[-1]['position']['phase'] == 'final_properties_on_saved_orbitals'
assert len(kept) == 3, [m['position'] for m in kept]
lengths = {rung(m)['index']: rung(m)['iterations_at_rung'] for m in ends}
iterations = result(fresh)['iterations_this_run']
assert sum(lengths.values()) == iterations, (lengths, iterations)
# The stop manifest carries the whole per-iteration history (A3, B2), one
# segment per rung, each counted from iteration 0.
history = latest(fresh)['metrics']['convergence_report']['iteration_history']
assert len(history) == iterations, (len(history), iterations)
for row in history:
    assert {'measurements', 'energy', 'kinetic_trace', 'explicit_accuracy', 'configured_progress',
            'engine_check_passed'} <= set(row), sorted(row)
    for name in ('density_change', 'bsh_residual', 'energy_change', 'charge_density_change',
                 'spin_density_change', 'delta_kinetic_energy'):
        assert name in row['measurements'], (name, sorted(row['measurements']))
samples = [(row['kinetic_trace']['sample']['threshold'], row['kinetic_trace']['sample']['iteration'])
           for row in history]
assert samples == [(1e-4, i) for i in range(lengths[1])] + [(1e-6, i) for i in range(lengths[2])], samples
assert 'earlier_runs_iteration_history' not in latest(fresh)['metrics']['convergence_report']
# The settings table shows the cadence, and the record says the same.
settings = result(fresh)['convergence_settings']['checkpoints']
assert settings['every_iterations'] is None and not settings['keep_all'], settings
assert settings['every_iterations_set_by'] == 'engine_default', settings
assert 'Checkpoints: at every rung end and stop' in log
assert result(fresh)['checkpoint_warnings'] == []

# 2. K = 3, everything kept: every third iteration at a rung, plus rung ends and the stop.
every3 = work / 'Every3KeepAll'
log = run('every3-keep-all', every3, with_keys('checkpoint.every_iterations 3', 'checkpoint.keep_all true'))
published = generations(every3)
mid = [(rung(m)['index'], rung(m)['iterations_at_rung']) for m in published if not is_decision_point(m)]
# A rung that stops on iteration n publishes its in-progress samples 1..n-1.
expected = [(r, i) for r in (1, 2) for i in range(1, lengths[r]) if i % 3 == 0]
assert mid == expected, (mid, expected)
assert sum(is_decision_point(m) for m in published) == 3, [m['position'] for m in published]
settings = result(every3)['convergence_settings']['checkpoints']
assert settings['every_iterations'] == 3 and settings['keep_all'], settings
assert settings['every_iterations_set_by'] == 'user' and settings['keep_all_set_by'] == 'user'
assert 'Checkpoints: every 3 iterations within a rung (user)' in log

# 3. K = 3 with retention: the mid-rung generations are replaced as the stage goes on.
every3_retained = work / 'Every3'
run('every3', every3_retained, with_keys('checkpoint.every_iterations 3'))
kept = generations(every3_retained)
assert all(is_decision_point(m) for m in kept) and len(kept) == 3, [m['position'] for m in kept]
assert not [d for d in (every3_retained / 'gs/checkpoints').iterdir() if d.name.startswith('.deleting-')]

# 4. --continue from a kept mid-rung generation (crash simulated right after it).
every1 = work / 'Every1KeepAll'
run('every1-keep-all', every1, with_keys('checkpoint.every_iterations 1', 'checkpoint.keep_all true'))
ordered = generations(every1)
assert len(ordered) == iterations + 1, len(ordered)  # every sample plus the final properties
chosen = next(m for m in ordered if rung(m)['index'] == 2 and not is_decision_point(m))
inside = work / 'InsideRung2'
shutil.copytree(every1, inside)
store = inside / 'gs/checkpoints'
for m in ordered[[g['generation'] for g in ordered].index(chosen['generation']) + 1:]:
    shutil.rmtree(store / m['generation'])
(store / 'latest').write_text(json.dumps({'generation': chosen['generation']}))
before = {m['generation'] for m in generations(inside)}
run('resume-inside-rung2', inside, resume=True)
record = result(inside)['continuation']
done = rung(chosen)['iterations_at_rung']
assert record['action'] == 'resumed' and record['rung']['index'] == 2, record
assert record['summary'] == 'resumed at rung 2 (thresh 1e-06, k 8), iteration %d of that rung' % (done + 1), record
assert result(inside)['status'] == 'completed'
after = generations(inside)
names = {m['generation'] for m in after}
# Default retention in the resumed run: the mid-rung generations, including
# the one it resumed from, are gone; rung ends and the new stop remain.
assert chosen['generation'] not in names, sorted(names)
assert all(is_decision_point(m) for m in after), [m['position'] for m in after]
assert [m['generation'] for m in after if m['generation'] in before] == \
       [m['generation'] for m in ordered if m['generation'] in before and is_decision_point(m)]
final = latest(inside)
assert rung(final)['iterations_at_rung'] == done + rung(final)['iterations_this_run'], rung(final)
report = final['metrics']['convergence_report']
earlier = report['earlier_runs_iteration_history']
assert len(earlier) == len(chosen['metrics']['convergence_report']['iteration_history']), len(earlier)
assert len(report['iteration_history']) == result(inside)['iterations_this_run']

# 5. A deletion that fails is a warning, the same on every rank; the stage
# finishes. Same crash as above; one older rung-1 generation cannot be emptied.
blocked = work / 'BlockedDeletion'
shutil.copytree(every1, blocked)
store = blocked / 'gs/checkpoints'
for m in ordered[[g['generation'] for g in ordered].index(chosen['generation']) + 1:]:
    shutil.rmtree(store / m['generation'])
(store / 'latest').write_text(json.dumps({'generation': chosen['generation']}))
victim = next(m for m in ordered if rung(m)['index'] == 1 and not is_decision_point(m))
(store / victim['generation']).chmod(stat.S_IRUSR | stat.S_IXUSR)  # its files cannot be removed
try:
    run('blocked-deletion', blocked, resume=True)
    warnings = result(blocked)['checkpoint_warnings']
    assert warnings and victim['generation'] in ' '.join(warnings), warnings
    assert result(blocked)['status'] == 'completed'
    # Whatever is left is under a .deleting- name, never listed as a generation.
    assert victim['generation'] not in {m['generation'] for m in generations(blocked)}
finally:
    for d in store.iterdir():
        if victim['generation'] in d.name:
            d.chmod(stat.S_IRWXU)
print('PASS')
