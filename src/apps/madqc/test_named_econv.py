#!/usr/bin/env python3
"""econv is the energy-change target (D8): only a user-written econv is a request,
judged on saved energy-change evidence (B6) and, in engine mode, reported only."""
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
work = Path(tempfile.mkdtemp(prefix='named-econv-', dir=Path.cwd()))
print('artifacts:', work, flush=True)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           OMPI_MCA_hwloc_base_binding_policy='none')
env.setdefault('MAD_NUM_THREADS', '10')  # dd convention; override from the environment
launch = [] if a.ranks == 1 else [a.mpiexec, '--bind-to', 'none', '-n', str(a.ranks)]
with_econv = (Path(a.source) / 'named_h2.in').read_text().split('molresponse optical')[0]
assert 'econv 1e-6' in with_econv
without_econv = with_econv.replace('  econv 1e-6\n', '')


def run(label, root, text, resume=False):
    inp = work / (label + '.in')
    inp.write_text(text)
    log = work / (label + '.log')
    with log.open('w') as out:
        r = subprocess.run(launch + [binary, str(inp), '--prefix=' + str(root)] +
                           (['--continue'] if resume else []), cwd=work, env=env,
                           stdout=out, stderr=subprocess.STDOUT, timeout=600)
    assert r.returncode == 0, log.read_text()[-6000:]
    print(label, 'exit', r.returncode, flush=True)


def result(root):
    return json.loads((root / 'gs/results.json').read_text())


def manifest(root):
    store = root / 'gs/checkpoints'
    generation = json.loads((store / 'latest').read_text())['generation']
    return json.loads((store / generation / 'manifest.json').read_text())


def copy(source, label):
    target = work / label
    shutil.copytree(source, target)
    return target


# A user-written econv is the energy-change target. It is reported, not
# enforced: the engine check still stops the run (B8 engine mode).
fresh = work / 'WithEconv'
run('with-econv', fresh, with_econv)
report = result(fresh)['convergence_report']
energy = report['measurements']['energy_change']
assert energy['target'] == 1e-6 and energy['set_by'] == 'user', energy
assert energy['relation'] == 'less_than_or_equal' and energy['active_gate'] is False, energy
assert energy['satisfied'] == (energy['value'] <= 1e-6), energy
saved = report['saved_state_validation']['measurements']['energy_change']
assert saved['target'] == 1e-6 and saved['availability'] == 'measured', saved
assert report['requested_accuracy']['criteria']['energy_change']['status'] in ('satisfied', 'unmet'), report['requested_accuracy']
assert manifest(fresh)['targets']['econv'] == 1e-6
assert result(fresh)['targets']['econv'] == 1e-6

# The engine default econv is not a target: no target, no verdict.
plain = work / 'WithoutEconv'
run('without-econv', plain, without_econv)
report = result(plain)['convergence_report']
energy = report['measurements']['energy_change']
assert energy['target'] is None and energy['satisfied'] is None, energy
assert energy['set_by'] == 'engine_default_not_a_target', energy
assert 'energy_change' not in report['requested_accuracy']['criteria']
assert 'econv' not in manifest(plain)['targets'] and 'econv' not in result(plain)['targets']
# econv was never part of the engine check, so the solve itself is unchanged.
assert result(plain)['iterations_this_run'] == result(fresh)['iterations_this_run']
assert abs(result(plain)['properties']['energy']['value'] - result(fresh)['properties']['energy']['value']) < 1e-12

# Reuse: the saved energy change is the evidence, not the recorded number.
saved_change = result(fresh)['convergence_report']['saved_state_validation']['measurements']['energy_change']['value']
loose = '%.1e' % max(1e-2, 10 * saved_change)
reuse = copy(fresh, 'LooseReuse')
run('loose-reuse', reuse, with_econv.replace('econv 1e-6', 'econv ' + loose), resume=True)
assert result(reuse)['work'] == 'loaded' and result(reuse)['iterations_this_run'] == 0
assert result(reuse)['convergence_report']['requested_accuracy']['criteria']['energy_change']['status'] == 'satisfied'

# A tighter econv than the saved evidence shows: in engine mode it can't be
# enforced, so the saved state is loaded and the request reported unmet
# (Adrian, 2026-10-21) rather than resumed only to stop on the engine check.
tight = copy(fresh, 'TightRequest')
run('tight-request', tight, with_econv.replace('econv 1e-6', 'econv 1e-14'), resume=True)
assert result(tight)['work'] == 'loaded' and result(tight)['iterations_this_run'] == 0
request = result(tight)['convergence_report']['requested_accuracy']
energy_request = request['criteria']['energy_change']
assert energy_request['status'] == 'unmet' and request['enforced'] is False, request
assert energy_request['outcome'] == 'unmet, not enforceable in engine mode', request
assert energy_request['value'] == saved_change, request
assert 'energy_change' in result(tight)['continuation']['unmet_not_enforceable']

# A run without econv can be continued with one; evidence still decides.
later = copy(plain, 'LaterRequest')
run('later-request', later, with_econv.replace('econv 1e-6', 'econv ' + loose), resume=True)
assert result(later)['work'] == 'loaded'

# Without econv, continuation is unchanged: loaded with no iterations.
again = copy(plain, 'PlainReuse')
run('plain-reuse', again, without_econv, resume=True)
assert result(again)['work'] == 'loaded' and result(again)['iterations_this_run'] == 0
print('PASS')
