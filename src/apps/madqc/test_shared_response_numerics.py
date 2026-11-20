#!/usr/bin/env python3
"""Compare shared results with the curated, unchanged scientific oracles.

The original check.json supplies every numerical tolerance. Historical format,
iteration and cost assertions remain in qctest; this integration check records
counts/timings without requiring the historical number of updates.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def lookup(value, keys):
    for key in keys:
        value = value[key]
    return value


def shared_value(key, reference, ground, properties):
    if key == ['tasks', 0, 'properties', 'energy']:
        return ground['total_energy']['value']
    if key[2:4] == ['metadata', 'excited_states'] and key[-1] == 'omega':
        return properties['excitation_energy']['components'][f'es_root_{key[-2]:04d}']
    if key[2:4] == ['properties', 'response_properties']:
        kind = key[4]
        old = lookup(reference, key[:7])
        if kind == 'alpha':
            row = next(r for r in properties['polarizability']['values']
                       if abs(r['frequency'] - old['omega']) < 1e-12)
            assert row['status'] == 'complete'
            if key[7] == 'omega':
                return row['frequency']
            directions = old.get('directions', 'xyz' if len(old['alpha']) == 3 else 'z')
            return row['components'][directions[key[8]] + directions[key[9]]]
        if kind in ('beta', 'raman'):
            name = 'hyperpolarizability' if kind == 'beta' else 'polarizability_gradient'
            rows = properties[name]['values']
            row = next(r for r in rows if abs(r['frequency_b'] - old['freq_b']) < 1e-12
                       and abs(r['frequency_c'] - old['freq_c']) < 1e-12
                       and (kind == 'beta' or r['coordinate'] == old['C']))
            assert row['status'] == 'complete'
            component = old['A'] + old['B'][-1] + (old['C'][-1] if kind == 'beta' else '')
            return row['components'][component]
        if kind == 'tpa':
            row = properties['two_photon_absorption']['values'][key[6]]
            assert row['status'] == 'complete'
            return row['observables'][key[7]]
    if key[2:4] == ['properties', 'vibrational_analysis']:
        values = properties['vibrational_analysis']['values']
        if key[4] == 'frequencies_cm1':
            # The conversion constant is part of the original calculation's record.
            factor = reference['tasks'][1]['properties']['vibrational_analysis']['au2invcm']
            return values['frequencies'][key[-1]] * factor
        if key[4] == 'reducedmass':
            return values['reduced_mass'][key[-1]]
    if key[2:4] == ['properties', 'raman_spectra']:
        return lookup(properties['raman_spectrum']['values'], key[5:])
    raise AssertionError(f'No shared scientific mapping for {key}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True)
    p.add_argument('--source', required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--work', help='Read an existing run from this script without solving again')
    a = p.parse_args()
    case = Path(a.source).resolve() / 'src/examples/qc' / a.case
    check = json.loads((case / 'check.json').read_text())
    reference = json.loads((case / 'reference' / f'{a.case}.calc_info.json').read_text())
    work = Path(a.work) if a.work else Path(tempfile.mkdtemp(prefix=f'shared-numerics-{a.case}-'))
    print('artifacts:', work, flush=True)
    if not a.work:
        for item in case.iterdir():
            if item.name in ('reference', 'check.json', 'README.md'):
                continue
            if item.is_dir():
                shutil.copytree(item, work / item.name)
            else:
                shutil.copy2(item, work / item.name)
        env = dict(os.environ, MAD_NUM_THREADS=str(a.threads), OPENBLAS_NUM_THREADS='1',
                   OMP_NUM_THREADS='1', OMPI_MCA_hwloc_base_binding_policy='none')
        with (work / 'run.log').open('w') as output:
            result = subprocess.run([str(Path(a.binary).resolve()), '--wf=response',
                                     '--prefix=run', f'{a.case}.in'], cwd=work, env=env,
                                    stdout=output, stderr=subprocess.STDOUT, timeout=14400)
        assert result.returncode == 0, f'calculation failed; see {work}/run.log'
    root = work / 'run'
    records = json.loads((root / 'results.json').read_text())['calculations']
    ground = next(c for c in records if c['calculation_id'] == 'ground_state')['properties']
    response = next(c for c in records if c['calculation_id'] == 'response')
    comparisons = []
    for entry in check['checks']:
        key = entry['key']
        if not ('rtol' in entry or entry.get('tol', 0) > 0):
            continue
        expected = lookup(reference, key)
        assert isinstance(expected, (float, int)) and not isinstance(expected, bool)
        actual = shared_value(key, reference, ground, response['properties'])
        tolerance = abs(expected) * entry['rtol'] if 'rtol' in entry else entry['tol']
        error = abs(actual - expected)
        comparisons.append(dict(reference_key=key, actual=actual, reference=expected,
                                tolerance=tolerance, absolute_error=error))
        print(key, actual, expected, 'error', error, 'limit', tolerance, flush=True)
    assert comparisons, 'no scientific comparisons selected'
    counts = []
    for latest in (root / 'checkpoints/response/tasks').glob('*/latest'):
        record = json.loads((latest.parent / json.loads(latest.read_text())['generation'] / 'manifest.json').read_text())
        task = record['task']
        assert task['outcome']['complete'], task
        if 'iteration' in task:
            counts.append({'task': latest.parent.name, 'attempt_id': record['attempt_id'],
                           'generation': record['generation'],
                           'updates': task['iteration']['updates_this_attempt']})
    evidence = {'case': a.case, 'reference': str(case / 'reference'), 'comparisons': comparisons,
                'counts': counts, 'counts_scope': 'latest checkpoint solve per task; reused checkpoints retain their original attempt',
                'timings': str(root / 'timings.json')}
    (work / 'numerical-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    assert all(math.isfinite(c['actual']) and c['absolute_error'] <= c['tolerance'] for c in comparisons), evidence
    print('shared scientific comparisons passed:', a.case, flush=True)


if __name__ == '__main__':
    main()
