#!/usr/bin/env python3
"""--continue resumes at the right rung and the iteration cap counts across runs
(DESIGN B5, B6): a finished rung is never redone, an interrupted rung resumes
with its iterations counted, and a final rung that used its cap is loaded and
reported unmet until the cap is raised."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='named-continuation-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
deck = (Path(a.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
assert 'protocol [1e-4, 1e-6]' in deck and 'maxiter 30' in deck


def run(label, root, text=deck, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=600)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


def result(root):
    return json.loads((root / 'gs/results.json').read_text())


def generations(root):
    """Published generations of the gs stage, oldest first."""
    store = root / 'gs/checkpoints'
    found = []
    for directory in store.iterdir():
        if directory.is_dir() and (directory / 'manifest.json').exists():
            found.append(json.loads((directory / 'manifest.json').read_text()))
    # g-<run start time>-<pid>-<sequence>: order by run, then by sequence.
    def order(m):
        parts = m['generation'].split('-')
        return int(parts[1]), int(parts[-1])
    return sorted(found, key=order)


def interrupted_at(source, label, choose):
    """A copy of `source` as a crash right after publishing the chosen
    generation leaves it: `latest` points there, nothing newer exists."""
    root = work / label
    shutil.copytree(source, root)
    store = root / 'gs/checkpoints'
    ordered = [m for m in generations(root) if m['generation'].startswith('g-')]
    chosen = next(m for m in ordered if choose(m))
    keep = True
    for m in ordered:
        if not keep:
            shutil.rmtree(store / m['generation'])
        if m['generation'] == chosen['generation']:
            keep = False
    (store / 'latest').write_text(json.dumps({'generation': chosen['generation']}))
    return root, chosen


def rung(m):
    return m['position']['rung']


# A fresh run records its rung in every generation.
fresh = work / 'Fresh'
run('fresh', fresh)
ordered = generations(fresh)
assert all('rung' in m['position'] for m in ordered), [m['position'] for m in ordered]
assert {rung(m)['index'] for m in ordered} == {1, 2}
ends = [m for m in ordered if m['position'].get('phase') == 'after_final_rotation']
assert [rung(m)['stop'] for m in ends] == ['converged', 'converged'], [rung(m) for m in ends]
fresh_energy = result(fresh)['properties']['energy']['value']

# Interrupted inside rung 2: --continue resumes rung 2 with its iteration count
# and never runs rung 1 again.
inside, chosen = interrupted_at(
    fresh, 'InsideRung2',
    lambda m: rung(m)['index'] == 2 and m['position']['phase'] == 'after_orbital_update')
done = rung(chosen)['iterations_at_rung']
run('resume-inside-rung2', inside, resume=True)
record = result(inside)['continuation']
assert record['action'] == 'resumed' and record['rung']['index'] == 2, record
assert record['iterations_done_at_rung'] == done, record
assert record['summary'] == 'resumed at rung 2 (thresh 1e-06, k 8), iteration %d of that rung' % (done + 1), record
new = [m for m in generations(inside) if m['generation'] not in {g['generation'] for g in ordered}]
assert new and all(rung(m)['index'] == 2 for m in new), [rung(m) for m in new]
assert rung(new[-1])['iterations_at_rung'] == done + result(inside)['iterations_this_run']
assert result(inside)['status'] == 'completed'
print('energy difference after mid-rung restart (acceleration history rebuilt):',
      result(inside)['properties']['energy']['value'] - fresh_energy, flush=True)

# Rung 1 finished: it is never redone; --continue starts rung 2 at iteration 1.
after_rung1, chosen = interrupted_at(
    fresh, 'AfterRung1',
    lambda m: rung(m)['index'] == 1 and m['position']['phase'] == 'after_final_rotation')
run('resume-after-rung1', after_rung1, resume=True)
record = result(after_rung1)['continuation']
assert record['rung']['index'] == 2 and record['iterations_done_at_rung'] == 0, record
new = [m for m in generations(after_rung1) if m['generation'] not in {g['generation'] for g in ordered}]
assert new and all(rung(m)['index'] == 2 for m in new), [rung(m) for m in new]

# Per-rung cap (maxiter [30, 2]): rung 2 stops at its cap, unconverged.
capped = work / 'Capped'
run('capped', capped, deck.replace('maxiter 30', 'maxiter [30, 2]'), success=False)
last = generations(capped)[-1]
assert rung(last)['index'] == 2 and rung(last)['iterations_at_rung'] == 2, rung(last)
assert rung(last)['stop'] == 'iteration_limit' and last['stop_reason'] != 'converged'
assert result(capped)['status'] == 'unconverged'

# Same cap on --continue: the cap is used up, so the state is loaded and
# reported unmet instead of iterating again.
same_cap = work / 'SameCap'
shutil.copytree(capped, same_cap)
run('same-cap', same_cap, deck.replace('maxiter 30', 'maxiter [30, 2]'), resume=True, success=False)
record = result(same_cap)['continuation']
assert record['action'] == 'loaded' and result(same_cap)['iterations_this_run'] == 0, record
assert record['unmet'] == 'unmet: iteration limit reached at rung 2', record

# A raised cap gives rung 2 more iterations; the count includes the first run.
raised = work / 'RaisedCap'
shutil.copytree(capped, raised)
run('raised-cap', raised, deck.replace('maxiter 30', 'maxiter [30, 30]'), resume=True)
record = result(raised)['continuation']
assert record['action'] == 'resumed' and record['rung']['index'] == 2, record
assert record['iterations_done_at_rung'] == 2, record
final = generations(raised)[-1]
assert rung(final)['iterations_at_rung'] == 2 + result(raised)['iterations_this_run'], rung(final)
assert result(raised)['status'] == 'completed'
print('PASS')
