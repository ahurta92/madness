#!/usr/bin/env python3
"""Who decides convergence (DESIGN B8) and the settings table (B9), for a named
moldft stage: the criteria-modes case table in both modes, the warnings, reuse
across modes, and that the printed table matches the record."""
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
work = named_workdir.work_directory('convergence-modes-')
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
base = (Path(a.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
assert '  dconv 1e-5\n' in base and '  econv 1e-6\n' in base
CRITERIA = '  convergence.use_default_criteria false\n'


def deck(*, criteria=False, dconv='1e-5', econv='1e-6', maxiter='30', extra=''):
    """The H2 moldft deck with the given written values (None: not written)."""
    text = base.replace('  dconv 1e-5\n', '' if dconv is None else '  dconv %s\n' % dconv)
    text = text.replace('  econv 1e-6\n', '' if econv is None else '  econv %s\n' % econv)
    text = re.sub(r'(?m)^\s*maxiter\s+\d+', '  maxiter ' + maxiter, text)
    return text.replace('moldft gs\n', 'moldft gs\n' + (CRITERIA if criteria else '') + extra)


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
                           stdout=out, stderr=subprocess.STDOUT, timeout=600)
    assert (r.returncode == 0) == success, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)
    stage = named_records.stages(root)[0]
    return root, stage, log.read_text()


def manifest(root):
    store = root / 'gs/checkpoints'
    return json.loads((store / json.loads((store / 'latest').read_text())['generation']
                       / 'manifest.json').read_text())


def table_rows(log):
    """The printed settings table, row name -> (target, set by, enforced)."""
    rows, inside = {}, False
    for line in log.splitlines():
        if line.startswith('Convergence settings:'):
            inside = True
        elif inside and line.startswith('  Iteration cap'):
            break
        elif inside:
            match = re.match(r'^    (.+?)\s{2,}(\S+)\s{2,}(engine default|user|not set)\s{2,}(yes|no)\s', line)
            if match:
                rows.setdefault(match.group(1), []).append(match.groups()[1:])
    return rows


def check_table_matches_record(log, stage):
    """B9 one source: every record row appears in the printed table as recorded."""
    rows = table_rows(log)
    settings = stage['convergence_report']['convergence_settings']
    assert ('Decided by: engine check' in log) == (settings['decided_by'] == 'engine_check')
    for metric in settings['metrics']:
        target = '-' if metric['target'] is None else '%.1e' % metric['target']
        expected = (target, metric['set_by'].replace('_', ' '), 'yes' if metric['enforced'] else 'no')
        assert expected in rows.get(metric['name'], []), (metric, rows)


# Row 1: both pass. The default deck (engine mode) completes as before.
root, stage, log = run('engine-default', deck())
assert stage['status'] == 'completed' and stage['convergence_report']['convergence_settings']['decided_by'] == 'engine_check'
check_table_matches_record(log, stage)
assert 'WARNING: dconv 1.0e-05 is reported, not enforced' in log
assert 'WARNING: econv 1.0e-06 is reported, not enforced' in log
assert stage['convergence_report']['convergence_settings']['warnings'] == [l.split('WARNING: ', 1)[1] for l in log.splitlines()
                                                     if l.startswith('  WARNING: ')]
engine_default = root

# A deck with no written targets at all shows only engine defaults.
_, plain, plain_log = run('engine-no-targets', deck(dconv=None, econv=None))
assert not [m for m in plain['convergence_report']['convergence_settings']['metrics']
            if m['group'] == 'requested_criteria' and m['target'] is not None]
assert 'engine default econv 1.0e-05 is not a target' in plain_log
check_table_matches_record(plain_log, plain)

# Row 1 in criteria mode: loose written targets, met.
_, stage, log = run('criteria-loose', deck(criteria=True, dconv='1e-3', econv='1e-3'))
assert stage['status'] == 'completed'
assert stage['convergence_report']['convergence_settings']['decided_by'] == 'requested_criteria'
check_table_matches_record(log, stage)
assert 'set no bsh_residual target' in log  # incomplete request warned (B8)

# Row 2: engine passes, written target unmet (econv 1e-14).
_, stage, log = run('engine-unmet-target', deck(econv='1e-14'))
assert stage['status'] == 'completed'
criterion = stage['convergence_report']['requested_accuracy']['criteria']['energy_change']
assert criterion['status'] == 'unmet' and criterion['outcome'] == 'unmet, not enforceable in engine mode'
root, stage, log = run('criteria-unmet-target', deck(criteria=True, econv='1e-14', maxiter='4'),
                       success=False)
assert stage['status'] == 'unconverged' and manifest(root)['stop_reason'] == 'iteration_limit'

# Row 3: engine fails, written target passes. econv 1 is met from the first
# comparable iteration; the engine check is not yet.
_, engine_run, _ = run('engine-loose-energy', deck(dconv=None, econv='1'))
root, criteria_run, log = run('criteria-loose-energy', deck(criteria=True, dconv=None, econv='1'))
assert criteria_run['status'] == 'completed'
assert criteria_run['iterations_this_run'] < engine_run['iterations_this_run'], \
    (criteria_run['iterations_this_run'], engine_run['iterations_this_run'])
history = manifest(root)['metrics']['convergence_report']['iteration_history']
assert history[-1]['engine_check_passed'] is False  # stopped by the criteria alone
assert manifest(root)['metrics']['convergence_report']['engine_check_passed'] is False

# Row 4: neither passes (one iteration per rung is never enough to compare).
_, stage, _ = run('engine-neither', deck(maxiter='1'), success=False)
assert stage['status'] == 'unconverged'
_, stage, _ = run('criteria-neither', deck(criteria=True, maxiter='1'), success=False)
assert stage['status'] == 'unconverged'

# Row 5: criteria mode with nothing written is rejected before solving.
inp = work / 'criteria-empty.in'
inp.write_text(deck(criteria=True, dconv=None, econv=None))
r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(work / 'CriteriaEmpty')], cwd=work,
                   env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
assert r.returncode != 0 and 'needs at least one written target' in r.stdout, r.stdout[-2000:]
assert not (work / 'CriteriaEmpty' / 'gs' / 'checkpoints').exists()

# Reuse across modes (B6, B8): an engine-mode generation satisfies a
# criteria-mode request only when its saved metrics meet every target.
saved = manifest(engine_default)['metrics']['convergence_report']['saved_state_validation']['measurements']
loose_energy = '%.1e' % max(1e-3, 10 * saved['energy_change']['value'])
_, stage, _ = run('criteria-reuse-met', deck(criteria=True, dconv=None, econv=loose_energy),
                  resume_from=engine_default)
assert stage['work'] == 'loaded' and stage['iterations_this_run'] == 0, stage['continuation']
spent = manifest(engine_default)['position']['protocol_step']['iterations_at_protocol_step']
_, stage, _ = run('criteria-reuse-unmet',
                  deck(criteria=True, dconv=None, econv='1e-14', maxiter=str(spent + 2)),
                  resume_from=engine_default, success=False)
assert stage['work'] == 'resumed' and stage['iterations_this_run'] > 0, stage['continuation']
print('PASS')
