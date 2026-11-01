#!/usr/bin/env python3
"""The stage handoff is disk-only (DESIGN A1): a molresponse stage reads its
reference only from the moldft stage's published generation (the orbital
archive and the consumer context in its manifest). Running gs alone, then a
deck holding only the response stage on that output directory, gives the same
results as running both stages in one deck."""
import argparse
import json
import sys
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
sys.dont_write_bytecode = True  # nothing written into the source tree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import named_records  # noqa: E402  (D-105: records by path)

p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--source', required=True)
p.add_argument('--ranks', type=int, default=1)
p.add_argument('--mpiexec', default='mpiexec')
a = p.parse_args()
binary = str(Path(a.binary).resolve())
work = Path(tempfile.mkdtemp(prefix='named-reference-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
# Off-centre and free to orient, so the SCF moves the molecule: the response
# must use the geometry the SCF used, which only the published context holds.
both = (Path(a.source) / 'named_h2.in').read_text()
assert '  no_orient true\n' in both and '  H 0 0 -0.7\n  H 0 0 0.7\n' in both
both = both.replace('  no_orient true\n', '').replace(
    '  H 0 0 -0.7\n  H 0 0 0.7\n', '  H 0.3 0 -0.4\n  H 0.3 0 1.0\n')
ground, response = both.split('molresponse optical')
response = 'molresponse optical' + response


def run(label, root, text, resume=False, success=True):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=900)
    if success is not None:  # None: either outcome is valid; the caller checks
        assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    return log.read_text()


def result(root, stage):
    return named_records.stage_record(root / stage)


def manifest(root, stage):
    store = root / stage / 'checkpoints'
    generation = json.loads((store / 'latest').read_text())['generation']
    return json.loads((store / generation / 'manifest.json').read_text())


def close(x, y, tol=1e-10):
    return abs(x - y) <= tol * max(1.0, abs(x), abs(y))


# One deck, both stages.
together = work / 'Together'
run('together', together, both)

# gs alone, then the response stage alone on the same output directory.
split = work / 'Split'
run('gs-alone', split, ground)
gs_generation = manifest(split, 'gs')['generation']
run('optical-alone', split, response, resume=True)

# The consumer context: the molecule the SCF used, not the deck's.
context = manifest(split, 'gs')['numerical']['consumer_context']
assert context['schema'] == 'madqc.moldft.consumer_context.v1', context
assert context['charge'] == 0 and context['multiplicity'] == 1, context
assert context['spin_restricted'] is True and context['functional'] == 'hf', context
assert context['thresh'] == 1e-6 and context['L'] == 20, context
geometry = result(split, 'gs')['geometry']['geometry']
assert context['molecule']['geometry'] == geometry, (context['molecule'], geometry)
# The SCF centred and oriented the molecule (bond along an axis, origin between).
assert all(abs(sum(atom[i] for atom in geometry)) < 1e-10 for i in range(3)), geometry
assert sorted(abs(c) for c in geometry[0]) [-1] == max(abs(c) for c in geometry[1]), geometry
assert any(abs(c) > 0.6 for c in geometry[0]) and geometry != [[0.3, 0, -0.4], [0.3, 0, 1.0]], geometry
assert context == manifest(together, 'gs')['numerical']['consumer_context']

# The response bound the published generation and gave the two-stage results.
alone, joint = result(split, 'optical'), result(together, 'optical')
assert alone['status'] == joint['status'] == 'completed', (alone['status'], joint['status'])
assert alone['reference'] == {'stage': 'gs', 'state_id': gs_generation}, alone['reference']
assert alone['reference_representation'] == joint['reference_representation']
assert alone['reference_representation'] == {'L': 20.0, 'k': context['k']}
assert alone['convergence_report']['reference']['state_id'] == gs_generation

def energy_evidence(record):
    """The reference's energy-change evidence, less what names the run."""
    evidence = dict(record['convergence_report']['reference']['energy_change'])
    evidence['sample'] = {k: v for k, v in evidence['sample'].items() if k != 'run_id'}
    return evidence.pop('value'), evidence


assert alone['convergence_report']['reference']['converged'] is True
assert joint['convergence_report']['reference']['converged'] is True
(value, evidence), (joint_value, joint_evidence) = energy_evidence(alone), energy_evidence(joint)
assert evidence == joint_evidence, (evidence, joint_evidence)
assert abs(value - joint_value) < 1e-12, (value, joint_value)  # gs ran twice: roundoff
for x, y in zip(alone['properties']['polarizability']['values'],
                joint['properties']['polarizability']['values'], strict=True):
    assert x['frequency'] == y['frequency'] and x['converged'] and y['converged']
    assert close(x['components']['zz'], y['components']['zz']), (x, y)


def iterations(record):
    return [(c['channel'], c['frequency'], c['protocol_step']['thresh'], c['iterations_at_save'])
            for c in record['convergence_report']['channels']]


assert iterations(alone) == iterations(joint), (iterations(alone), iterations(joint))

# A reference outside the deck is read only from a published generation; one
# with none is refused before any stage runs and before a fresh root is
# created (D-100).
log = run('optical-fresh', work / 'NoReference', response, success=False)
assert 'reference gs names no stage of this deck (optical) and no published generation' in log, log[-3000:]
assert not (work / 'NoReference').exists()
# Tightening is checked against the published generation, too.
tight = response.replace('protocol [1e-4, 1e-6]', 'protocol [1e-4, 1e-8]')
log = run('optical-tight', split, tight, resume=True, success=False)
# D-102: the refusal names both sides and asks for the protocol meant.
assert "response stage optical: protocol final step thresh 1.0e-08" in log, log[-3000:]
assert "is tighter than reference gs's saved generation g-gs-" in log, log[-3000:]
assert 'write the protocol you mean' in log

# A reused ground state is the same reference in and out of the deck (A1, B6).
# Run 1: the written targets decide and cannot be met in 3 iterations per
# rung, so the final rung stops at its iteration limit.
assert '  econv 1e-6\n' in ground and '  dconv 1e-5\n' in ground and '  maxiter 30\n' in ground
capped = ground.replace('moldft gs\n', 'moldft gs\n  convergence.use_default_criteria false\n') \
               .replace('  econv 1e-6\n', '  econv 1e-14\n').replace('  maxiter 30\n', '  maxiter 3\n')
reused = work / 'Reused'
run('gs-capped', reused, capped, success=False)
before = manifest(reused, 'gs')
assert before['stop_reason'] == 'iteration_limit', before['stop_reason']
# Run 2: looser targets that the saved orbitals meet. The state is loaded, and
# because its generation did not stop converged, the verdict is published as a
# new generation with the same orbital shards.
loose = capped.replace('  econv 1e-14\n', '  econv 1e-2\n').replace('  dconv 1e-5\n', '  dconv 1e-1\n')
run('gs-reused', reused, loose, resume=True)
record = result(reused, 'gs')
assert record['work'] == 'loaded' and record['iterations_this_run'] == 0, record['work']
assert record['status'] == 'completed', record['status']
after = manifest(reused, 'gs')
assert after['generation'] != before['generation'] and record['state_id'] == after['generation']
assert after['stop_reason'] == 'converged' and after['position']['reused_from'] == before['generation']
verdict = after['metrics']['convergence_report']['reuse_verdict']
assert verdict['stop_reason_before'] == 'iteration_limit' and verdict['decided_by'] == 'requested_criteria', verdict
assert [s['checksum'] for s in after['archives'][0]['shards']] == \
       [s['checksum'] for s in before['archives'][0]['shards']]
assert after['position']['protocol_step'] == before['position']['protocol_step']  # how the rung ended stays recorded
# Run 3: a response stage outside gs's deck now accepts that generation, as an
# in-deck reference does.
run('optical-on-reused', reused, response, resume=True, success=None)
assert result(reused, 'optical')['reference']['state_id'] == after['generation']
print('disk-only handoff: response from the published generation equals the two-stage run')
