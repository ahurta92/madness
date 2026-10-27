#!/usr/bin/env python3
"""Exercise response state reuse independently from result-file presence."""
import argparse
import hashlib
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
p.add_argument('--identity-only', action='store_true')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='response-checkpoint-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '-n', str(a.ranks)]
deck = (Path(a.source) / 'named_h2.in').read_text()
one = deck.replace('[0.0, 0.02]', '[0.0]')

def run(label, root, text=one, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=300)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()

def published(root, stage='optical'):
    store = root / stage / 'checkpoints'
    assert (store / 'latest').exists(), 'response did not publish a checkpoint'
    directory = store / json.loads((store / 'latest').read_text())['generation']
    return directory, json.loads((directory / 'manifest.json').read_text())

def result(root):
    return json.loads((root / 'optical/results.json').read_text())

def hashes(directory):
    return {str(f.relative_to(directory)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in directory.rglob('*') if f.is_file()}

def alpha(root, index=0):
    return result(root)['properties']['polarizability']['values'][index]['components']['zz']

def check_identity(base):
    failures = []
    # Engine filenames round these different values to the same keys.
    for label, changed in (
            ('frequency-alias', deck.replace('0.02]', '0.020001]')),
            ('protocol-alias', '[1e-4, 1.01e-6]'.join(deck.rsplit('[1e-4, 1e-6]', 1)))):
        root = work / label
        shutil.copytree(base, root)
        run(label, root, changed, resume=True)
        if result(root)['iterations_this_run'] == 0:
            failures.append(label + ': changed scientific state was reused')
    root = work / 'IndependentSubset'
    shutil.copytree(base, root)
    tight_subset = 'dconv 1e-5'.join(one.rsplit('dconv 1e-4', 1))
    run('tight-subset', root, tight_subset, resume=True)
    independent = deck.replace('[0.0, 0.02]', '[0.02]')
    run('independent-subset', root, independent, resume=True)
    if result(root)['iterations_this_run'] != 0:
        failures.append('tightening a subset discarded independent convergence evidence')
    assert not failures, '; '.join(failures)

if a.identity_only:
    base = work / 'IdentityBase'
    run('identity-base', base, deck)
    check_identity(base)
    print('exact identity and independent reuse passed', flush=True)
    raise SystemExit(0)

fresh = work / 'Fresh'
run('fresh', fresh)
generation, manifest = published(fresh)
assert manifest['stop_reason'] == 'converged'
assert manifest['archives']
assert manifest['reference']['generation'] == published(fresh, 'gs')[1]['generation']
assert all(p['generation'] == manifest['generation'] for p in manifest['properties'])
assert result(fresh)['iterations_this_run'] > 0
fresh_result = result(fresh)
original = hashes(generation)
value = alpha(fresh)
assert abs(value - 6.4495330229835455) < 2e-5
# Results are derived output; checkpoint reuse must survive their removal.
(fresh / 'optical/results.json').unlink()
run('reuse', fresh, resume=True)
assert result(fresh)['work'] == 'loaded'
assert result(fresh)['iterations_this_run'] == 0
assert published(fresh)[0] == generation and hashes(generation) == original
assert alpha(fresh) == value

baseline = work / 'Baseline'
shutil.copytree(fresh, baseline)

# An additional frequency preserves the completed independent state archives.
run('add-frequency', fresh, deck, resume=True)
added, added_manifest = published(fresh)
assert result(fresh)['iterations_this_run'] > 0
assert abs(alpha(fresh) - value) < 1e-12
assert abs(alpha(fresh, 1) - 6.459428018020927) < 2e-5
for path, digest in original.items():
    if path != 'manifest.json':
        assert hashes(added)[path] == digest
assert hashes(generation) == original
# R4/R7: the new frequency starts at rung 1 from the nearest saved one (omega 0).
added_states = {s['frequency']: s for s in result(fresh)['continuation']['states']}
assert added_states[0.0]['action'] == 'loaded'
assert added_states[0.02]['action'] == 'seeded' and added_states[0.02]['seed_frequency'] == 0.0
assert added_states[0.02]['rung']['index'] == 1
assert [i['work'] for i in result(fresh)['work_items']][0] == 'seeded'
print('add-frequency iterations (seeded from omega 0):',
      [(i['protocol'], i['iterations']) for i in result(fresh)['work_items']], flush=True)

# R5: a frequency dropped from the request is not planned, loaded or deleted.
removed = work / 'Removed'
shutil.copytree(fresh, removed)
run('removed-frequency', removed, one, resume=True)
assert [s['frequency'] for s in result(removed)['continuation']['states']] == [0.0]
assert result(removed)['iterations_this_run'] == 0
kept = hashes(published(removed)[0])
assert any('f0.02000' in path for path in kept)
for path, digest in hashes(added).items():
    if 'f0.02000' in path:
        assert kept[path] == digest, path

# The planner decides each state on its own (DESIGN B6, D-016). Selecting the
# generation saved after omega 0.02 finished rung 1 leaves omega 0 complete at
# the final rung: it loads, and omega 0.02 resumes at rung 2 only.
per_state = work / 'PerState'
shutil.copytree(fresh, per_state)
store = per_state / 'optical/checkpoints'
selected = [m for m in (json.loads(f.read_text()) for f in store.glob('g-*/manifest.json'))
            if m['stop_reason'] == 'in_progress' and m['position']['frequency'] == 0.02]
selected = min(selected, key=lambda m: int(m['generation'].rsplit('-', 1)[1]))
first_rung = selected['position']['protocol']
(store / 'latest').write_text(json.dumps({'generation': selected['generation']}))
run('per-state', per_state, deck, resume=True)
states = {s['frequency']: s for s in result(per_state)['continuation']['states']}
assert states[0.0]['action'] == 'loaded' and states[0.0]['rung']['index'] == 2
assert states[0.02]['action'] == 'resumed' and states[0.02]['rung']['index'] == 2
items = result(per_state)['work_items']
assert [(i['frequency'], i['protocol'] != first_rung) for i in items] == [(0.02, True)]
assert abs(alpha(per_state) - alpha(fresh)) < 1e-12
assert abs(alpha(per_state, 1) - alpha(fresh, 1)) < 1e-8
# D-016: the loaded leg is met by the planner's verdict, and its row says so;
# the engine's own flag on it is unchanged.
values = result(per_state)['properties']['polarizability']['values']
loaded_row = values[0]['engine_evidence']
assert values[0]['converged'] and loaded_row['decided_by'] == 'planner'
assert all(leg['met'] and leg['decided_by'] == 'planner' and leg['converged']
           for leg in loaded_row['row_accuracy'])
assert 'decided_by' not in values[1]['engine_evidence']
# The adapter no longer edits the engine's file (D7): every state this run did
# not solve is in the attempt exactly as the selected generation saved it.
attempt = Path(published(per_state)[1]['numerical']['attempt_directory'])
after = json.loads((attempt / 'response_metadata.json').read_text())['fd_states']
saved = selected['numerical']['engine_metadata']['fd_states']
for channel, rungs in saved.items():
    for key, points in rungs.items():
        for fkey, entry in points.items():
            if not (fkey == 'f0.02000' and key != first_rung):
                assert after[channel][key][fkey] == entry, (channel, key, fkey)

partial = work / 'Partial'
short = one.rsplit('maxiter 30', 1)
short = 'maxiter 1'.join(short)
run('partial', partial, short, success=False)
old_dir, old = published(partial)
assert old['stop_reason'] != 'converged'
assert result(partial)['status'] == 'unconverged'
old_hash = hashes(old_dir)
run('resume', partial, resume=True)
assert result(partial)['work'] == 'resumed'
assert result(partial)['iterations_this_run'] > 0
assert abs(alpha(partial) - value) < 2e-4
assert hashes(old_dir) == old_hash

# All required members must be checked before entering collective readers.
for defect in ('missing', 'corrupt', 'inventory', 'definition', 'prerequisite'):
    broken = work / defect
    shutil.copytree(fresh, broken)
    directory, m = published(broken)
    shard = directory / m['archives'][0]['shards'][0]['path']
    if defect == 'missing':
        shard.unlink()
    elif defect == 'corrupt':
        data = bytearray(shard.read_bytes()); data[-1] ^= 1; shard.write_bytes(data)
    else:
        states = m['numerical']['engine_metadata']['fd_states']
        if defect == 'inventory':
            m['archives'].pop()
        elif defect == 'definition':
            entry = next(iter(next(iter(next(iter(states.values())).values())).values()))
            entry['type'] = 'unsupported-model'
        else:
            removed = set()
            for channel in states.values():
                for rung in channel.values():
                    for key in list(rung):
                        if rung[key]['freq'] == 0.0:
                            removed.add(rung.pop(key)['archive'])
            m['archives'] = [entry for entry in m['archives'] if entry['basename'] not in removed]
        (directory / 'manifest.json').write_text(json.dumps(m))
    before = hashes(broken / 'optical/checkpoints')
    run(defect, broken, deck, resume=True, success=False)
    assert hashes(broken / 'optical/checkpoints') == before

# A complete state bundle interrupted before analysis needs only analysis.
analysis = work / 'Analysis'
shutil.copytree(baseline, analysis)
store = analysis / 'optical/checkpoints'
intermediate = sorted((p for p in store.glob('g-*/manifest.json')
                       if json.loads(p.read_text())['stop_reason'] == 'in_progress'),
                      key=lambda p: int(p.parent.name.rsplit('-', 1)[1]))[-1]
(store / 'latest').write_text(json.dumps({'generation': intermediate.parent.name}))
run('analysis-only', analysis, resume=True)
assert result(analysis)['work'] == 'analyzed'
assert result(analysis)['iterations_this_run'] == 0
assert abs(alpha(analysis) - value) < 1e-12

# Unpublished files and failed later publication cannot replace the prior state.
lock_root = work / 'PublicationFailure'
shutil.copytree(baseline, lock_root)
old_dir, old = published(lock_root)
old_hash = hashes(old_dir)
(lock_root / 'optical/checkpoints/unpublished').mkdir()
(lock_root / 'optical/checkpoints/unpublished/junk').write_bytes(b'incomplete')
run('ignore-unpublished', lock_root, resume=True)
lock = lock_root / 'optical/checkpoints/.publish-lock'
lock.mkdir()
run('failed-publication', lock_root, deck, resume=True, success=False)
assert published(lock_root)[0] == old_dir and hashes(old_dir) == old_hash
lock.rmdir()
run('reuse-after-failure', lock_root, resume=True)
assert result(lock_root)['iterations_this_run'] == 0

# Tighter response targets must not inherit a looser convergence claim.
tighter = work / 'Tighter'
shutil.copytree(baseline, tighter)
tight_deck = one.rsplit('dconv 1e-4', 1)
run('tighter-response', tighter, 'dconv 1e-5'.join(tight_deck), resume=True)
assert result(tighter)['work'] == 'resumed' and result(tighter)['iterations_this_run'] > 0
assert abs(alpha(tighter) - value) < 2e-3

# A coarser working response preserves the original fine ground archive.
coarse = work / 'CoarserResponse'
shutil.copytree(baseline, coarse)
gs_dir, _ = published(coarse, 'gs')
gs_hash = hashes(gs_dir)
coarse_deck = '[1e-4]'.join(one.rsplit('[1e-4, 1e-6]', 1))
run('coarser-response', coarse, coarse_deck, resume=True)
assert hashes(gs_dir) == gs_hash
unsupported = '[1e-4, 1e-7]'.join(one.rsplit('[1e-4, 1e-6]', 1))
run('unsupported-tightening', coarse, unsupported, resume=True, success=False)
assert hashes(gs_dir) == gs_hash

# Explicit upstream refinement produces a new generation and invalidates descendants.
changed = work / 'ChangedReference'
shutil.copytree(baseline, changed)
old_dir, old = published(changed)
old_hash = hashes(old_dir)
refined = one.replace('[1e-4, 1e-6]', '[1e-4, 1e-6, 1e-7]')
run('changed-reference', changed, refined, resume=True)
assert result(changed)['reference_changed']
assert result(changed)['work'] == 'solved'
assert result(changed)['iterations_this_run'] > 0
assert published(changed)[1]['reference'] != old['reference']
assert hashes(old_dir) == old_hash
assert abs(alpha(changed) - value) < 2e-4

check_identity(fresh)

# Newly added independent stages must have distinct state identities even when
# their channel/protocol counts are identical in the same invocation.
twins = work / 'IndependentStages'
shutil.copytree(baseline, twins)
ground, response_text = one.split('molresponse optical')
response_text = response_text.replace('[1e-4, 1e-6]', '[1e-4]')
run('independent-stages', twins, ground + 'molresponse first' + response_text +
    '\nmolresponse second' + response_text, resume=True)
rows = json.loads((twins / 'results.json').read_text())['stages']
assert len({row['state_id'] for row in rows}) == len(rows)

summary = {name: {'iterations': result(root)['iterations_this_run'],
                  'work': result(root)['work'], 'alpha_zz': alpha(root)}
           for name, root in [('reused', baseline), ('added_frequency', fresh),
                              ('resumed', partial), ('analysis', analysis),
                              ('tighter', tighter), ('changed_reference', changed)]}
summary['fresh'] = {'iterations': fresh_result['iterations_this_run'], 'work': fresh_result['work'],
                    'alpha_zz': value}
(work / 'summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
print('response checkpoint tests passed', flush=True)
