#!/usr/bin/env python3
"""Who decides convergence (DESIGN B8) and the settings table (B9) for a named
molresponse stage (H2, static alpha): engine check vs written targets, the
warnings, a criteria request without targets, reuse across modes, and the
printed table against the record."""
import argparse
import json
import sys
import os
from pathlib import Path
import re
import shutil
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
work = named_workdir.work_directory('response-modes-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
# One rung and one frequency: these cases test who decides, not the ladder.
base = (Path(a.source) / 'named_h2_plumbing.in').read_text()
assert base.count('  dconv 1e-4\n') == 1  # the response stage's dconv


def deck(*, criteria=False, dconv='1e-4', maxiter='30', extra=''):
    """H2 deck; the response stage gets the given written values."""
    head, response = base.split('molresponse optical\n', 1)
    response = response.replace('  dconv 1e-4\n', '' if dconv is None else '  dconv %s\n' % dconv)
    response = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter ' + maxiter, response)
    lines = ('  convergence.use_default_criteria false\n' if criteria else '') + extra
    return head + 'molresponse optical\n' + lines + response


def run(label, text, success=True, resume_from=None):
    root = work / label
    if resume_from:
        shutil.copytree(resume_from, root)
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume_from else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=900)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    stages = named_records.stages(root)
    return root, stages[-1], log.read_text()


def entries(root):
    store = root / 'optical/checkpoints'
    m = json.loads((store / json.loads((store / 'latest').read_text())['generation']
                    / 'manifest.json').read_text())
    for channel in m['numerical']['engine_metadata']['fd_states'].values():
        for rung in channel.values():
            yield from rung.values()


def check_table_matches_record(log, stage):
    """B9 one source: each recorded row is printed with its target, setter and enforcement."""
    section = log.split('Convergence settings: stage optical', 1)[1].split('Iteration cap', 1)[0]
    settings = stage['convergence_report']['convergence_settings']
    assert ('Decided by: engine check' in section) == (settings['decided_by'] == 'engine_check')
    for metric in settings['metrics']:
        target = '-' if metric['target'] is None else '%.1e' % metric['target']
        row = r'\n    %s\s{2,}%s\s{2,}%s\s{2,}%s\s' % (re.escape(metric['name']), re.escape(target),
                                                    metric['set_by'].replace('_', ' '),
                                                    'yes' if metric['enforced'] else 'no')
        assert re.search(row, section), (metric, section)


# Engine decides (default): as before; the written dconv is reported.
engine_root, stage, log = run('engine-default', deck())
assert stage['status'] == 'completed'
assert stage['convergence_report']['convergence_settings']['decided_by'] == 'engine_check'
assert stage['convergence_report']['decided_by'] == 'engine_check'
check_table_matches_record(log, stage)
assert 'WARNING: dconv 1.0e-04 is reported, not enforced' in log
assert 'Convergence outcome: stage optical (molresponse)' in log
engine_iterations = stage['iterations_this_run']

# Engine decides; a written target it doesn't enforce stays unmet and reported.
_, stage, _ = run('engine-unmet-target', deck(extra='  accuracy.kinetic_trace_change 0\n'))
assert stage['status'] == 'completed'
stage_root = work / 'engine-unmet-target'
assert any(e['saved_state_validation']['explicit_accuracy']['status'] == 'unmet'
           for e in entries(stage_root) if 'explicit_accuracy' in e.get('saved_state_validation', {}))
assert all(e.get('checkpoint_decided_by') == 'engine_check' for e in entries(stage_root))

# Targets decide, loose: met, decided by the targets.
_, stage, log = run('criteria-loose', deck(criteria=True, dconv='1e-2',
                                           extra='  accuracy.bsh_residual 1e-1\n'))
assert stage['status'] == 'completed'
assert stage['convergence_report']['convergence_settings']['decided_by'] == 'requested_criteria'
check_table_matches_record(log, stage)
assert all(e.get('checkpoint_decided_by') == 'requested_criteria' for e in entries(work / 'criteria-loose'))

# Targets decide, unreachable (kinetic 0): iterates to the cap, unconverged,
# no published properties.
root, stage, _ = run('criteria-unmet', deck(criteria=True, maxiter='3',
                                            extra='  accuracy.kinetic_trace_change 0\n'),
                     success=False)
assert stage['status'] == 'unconverged'
assert all(not e['converged'] and not e.get('accepted', False) for e in entries(root))

# Targets decide with nothing written: rejected before solving.
inp = work / 'criteria-empty.in'
inp.write_text(deck(criteria=True, dconv=None))
r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(work / 'CriteriaEmpty')], cwd=work,
                   env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
assert r.returncode != 0 and 'needs at least one written target' in r.stdout, r.stdout[-2000:]

# Reuse across modes (B8): the engine-mode response meets loose written targets
# from its saved evidence, so a criteria request needs no new iterations.
_, stage, _ = run('criteria-reuse', deck(criteria=True, dconv='1e-2'), resume_from=engine_root)
assert stage['iterations_this_run'] == 0, stage.get('work_items')
print('engine iterations', engine_iterations)
print('PASS')
