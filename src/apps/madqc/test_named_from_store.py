#!/usr/bin/env python3
"""The from-store tier, level 2 (D-079, D-080): continue stored generations
with --continue, so the read, the plan, the publish, retention and the
handoff to a later stage run from a load, with no solve.

The store is make_test_store.py's (test_store.json lists it). Without
MADQC_TEST_STORE the test skips with one line; with it, the store is
verified first (a mismatch fails, nothing re-solves), and a store written by
a commit that is not an ancestor of this checkout fails with the instruction
to regenerate it (D-084).

Cases (each copies a stored run into a scratch directory first):
1. the same deck: every stage loaded, 0 iterations, the stored verdict and
   properties, no new generation, no archive shards left in the attempt;
2. a response-only deck naming the stored ground state: the handoff from a
   published generation (A1);
3. the capped ground state with looser targets: the reuse republish
   (decided_by planner, reused_from, the shards byte for byte), one
   generation left (D-070);
4. a 1-rank run continued on 2 ranks and the reverse (FD states load across
   rank counts);
5. a frequency removed from the request (R5): nothing solved, the removed
   frequency's state kept.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()

if not os.environ.get('MADQC_TEST_STORE'):
    print('skipped: MADQC_TEST_STORE not set')
    sys.exit(77)

sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, a.source)
import make_test_store as store  # noqa: E402

data = store.listing()
ancestor = subprocess.run(['git', '-C', a.source, 'merge-base', '--is-ancestor',
                           data['commit'], 'HEAD'], capture_output=True)
assert ancestor.returncode == 0, (
    f"the store was written by {data['commit']}, which is not an ancestor of this "
    "checkout: regenerate it (make_test_store.py generate <madqc>)")
problems = store.verify()
assert not problems, 'the store does not match test_store.json:\n' + '\n'.join(problems[:20])
stored = store.store_dir(data)

env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')
work = Path(tempfile.mkdtemp(prefix='from-store-', dir=Path.cwd()))
print('artifacts:', work, flush=True)


def copy(name, label):
    root = work / label
    shutil.copytree(stored / name, root)
    return root


def run(label, root, text, ranks=1):
    inp = work / (label + '.in')
    inp.write_text(text)
    launch = [] if ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(ranks)]
    with (work / (label + '.log')).open('w') as out:
        r = subprocess.run(launch + [a.binary, str(inp), f'--prefix={root}', '--continue'],
                           cwd=work, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=600)
    assert r.returncode == 0, (work / (label + '.log')).read_text()[-4000:]
    print(label, 'exit 0', flush=True)


def result(root, stage):
    return json.loads((root / stage / 'results.json').read_text())


def latest(root, stage):
    store_dir = root / stage / 'checkpoints'
    name = json.loads((store_dir / 'latest').read_text())['generation']
    return name, json.loads((store_dir / name / 'manifest.json').read_text())


def generations(root, stage):
    return sorted(d.name for d in (root / stage / 'checkpoints').iterdir()
                  if d.is_dir() and not d.name.startswith('.'))


def no_shards(root, stage, run_id):
    attempt = root / stage / ('attempt-' + run_id)
    shards = [f for f in attempt.rglob('*') if f.is_file() and f.suffix[1:].isdigit()]
    assert not shards, (stage, shards)


def deck(run_name):
    return (stored / f'{run_name}.in').read_text()


def run_id(root):
    return json.loads((root / 'results.json').read_text())['run_id']


def loaded_as_stored(root, run_name, stages):
    for stage in stages:
        before, after = result(stored / run_name, stage), result(root, stage)
        assert after['iterations_this_run'] == 0, (stage, after['iterations_this_run'])
        assert after['status'] == before['status'], (stage, after['status'])
        assert after['properties'] == before['properties'], stage
        assert latest(root, stage)[0] == latest(stored / run_name, stage)[0], stage
        no_shards(root, stage, run_id(root))


# 1. The same deck, at the ranks each run was written with.
for name, stages, ranks in (('h2-r1', ('gs', 'optical'), 1), ('he-r1', ('gs', 'resp'), 1),
                            ('h2-r2', ('gs', 'optical'), 2), ('he-r2', ('gs', 'resp'), 2)):
    root = copy(name, f'same-{name}')
    run(f'same-{name}', root, deck(name), ranks)
    assert result(root, 'gs')['work'] == 'loaded'
    loaded_as_stored(root, name, stages)

# 2. A response-only deck: the reference read from the stored ground state.
root = copy('h2-r1', 'response-only')
text = deck('h2-r1')
response_only = text.split('moldft gs')[0] + 'molresponse optical' + text.split('molresponse optical')[1]
run('response-only', root, response_only)
assert not (root / 'gs' / f"attempt-{run_id(root)}").exists(), 'the ground state ran'
loaded_as_stored(root, 'h2-r1', ('optical',))
assert latest(root, 'optical')[1]['reference']['generation'] == latest(stored / 'h2-r1', 'gs')[0]

# 3. The capped ground state, judged against looser targets: the verdict is
# published again with the same shards; retention leaves that one generation.
root = copy('h2-capped-r1', 'reused')
source, before = latest(root, 'gs')
assert before['stop_reason'] == 'iteration_limit', before['stop_reason']
loose = deck('h2-capped-r1').replace('  econv 1e-14\n', '  econv 1e-2\n').replace('  dconv 1e-5\n', '  dconv 1e-1\n')
run('reused', root, loose)
record = result(root, 'gs')
assert record['work'] == 'loaded' and record['iterations_this_run'] == 0, record['work']
assert record['status'] == 'completed', record['status']
name, after = latest(root, 'gs')
assert name != source and after['position']['reused_from'] == source
assert after['stop_reason'] == 'converged'
assert after['numerical']['states'][0]['decided_by'] == 'planner'
assert [s['checksum'] for s in after['archives'][0]['shards']] == \
       [s['checksum'] for s in before['archives'][0]['shards']]
assert generations(root, 'gs') == [name], generations(root, 'gs')
no_shards(root, 'gs', run_id(root))

# 4. Across rank counts: a 1-rank store continued on 2 ranks, and the reverse.
for name, ranks in (('h2-r1', 2), ('h2-r2', 1)):
    root = copy(name, f'ranks-{name}-on-{ranks}')
    run(f'ranks-{name}-on-{ranks}', root, deck(name), ranks)
    loaded_as_stored(root, name, ('gs', 'optical'))

# 5. A frequency removed from the request (R5): planned no more, solved no
# more, its state kept in what the stage publishes.
root = copy('he-r1', 'removed-frequency')
run('removed-frequency', root, deck('he-r1').replace('[0.0, 0.05]', '[0.0]'))
assert result(root, 'resp')['iterations_this_run'] == 0
assert [s['frequency'] for s in result(root, 'resp')['continuation']['states']] == [0.0]
name, manifest = latest(root, 'resp')
assert any('f0.05000' in s['key'] for s in manifest['numerical']['states']), \
    [s['key'] for s in manifest['numerical']['states']]
assert len(generations(root, 'resp')) == 1, generations(root, 'resp')

print('from-store tests passed', flush=True)
