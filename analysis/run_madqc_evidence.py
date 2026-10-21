#!/usr/bin/env python3
"""Plan and run isolated MADQC convergence evidence cases."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

SOURCE = Path('/home/ahurtado/repos/madness-step4-astra')
DEFAULT_MADQC = Path('/tmp/madness-hdf5-analysis/src/apps/madqc/madqc')
DEFAULT_OUT = Path('/tmp/madqc-evidence-runs')

@dataclass(frozen=True)
class Case:
    case_id: str
    molecule: str
    deck: Path
    workflow: str
    response_kind: str
    protocols: tuple[float, ...]
    forced_budgets: tuple[int, ...]

CASES = (
    Case('h2', 'H2', SOURCE/'src/examples/qc/response_h2_es_rpa_nokain/response_h2_es_rpa_nokain.in', 'response', 'static+excited', (1e-4, 1e-6), (5, 10)),
    Case('lih', 'LiH', SOURCE/'src/examples/qc/response_lih_beta/response_lih_beta.in', 'response', 'static', (1e-4,), (5, 10)),
    Case('water', 'H2O', SOURCE/'src/examples/qc/response_h2o_lda_alpha/response_h2o_lda_alpha.in', 'response', 'static', (1e-4, 1e-6), (5, 10)),
)


def plan_records(out_root: Path) -> list[dict]:
    records = []
    for case in CASES:
        for threshold in case.protocols:
            endpoint = f'{case.case_id}-gs-{threshold:.0e}'
            for budget in case.forced_budgets:
                records.append({
                    'run_id': f'{endpoint}-forced-{budget}', 'molecule': case.molecule,
                    'workflow': 'moldft', 'source_branch': 'madqc/step4-astra',
                    'ground_state_endpoint': None,
                    'protocol': {'threshold': threshold, 'forced_budget': budget},
                    'phase': 'ground_state_endpoint', 'sampling_phase': 'exhausted_rung',
                    'classification': 'forced', 'artifacts': {'output_dir': str(out_root/endpoint/f'forced-{budget}')},
                    'status': 'planned',
                })
            records.append({
                'run_id': f'{case.case_id}-ordinary', 'molecule': case.molecule,
                'workflow': 'response', 'source_branch': 'madqc/step4-astra',
                'ground_state_endpoint': endpoint,
                'protocol': {'thresholds': list(case.protocols)},
                'phase': 'response_iteration', 'sampling_phase': 'ordinary',
                'classification': 'ordinary', 'artifacts': {'output_dir': str(out_root/case.case_id/'ordinary')},
                'status': 'planned',
            })
    return records


def prepare_deck(source: Path, destination: Path) -> None:
    text = source.read_text()
    if re.search(r'\bhdf5\s+false\b', text):
        text = re.sub(r'(\bhdf5\s+)false\b', r'\1true', text, count=1)
    elif not re.search(r'^\s*hdf5\s+true\b', text, flags=re.MULTILINE):
        text = re.sub(r'(^response\s*$)', r'\1\n    hdf5 true', text, count=1, flags=re.MULTILINE)
    destination.write_text(text)


def execute_case(case: Case, madqc: Path, out_root: Path) -> dict:
    case_out = out_root / case.case_id / 'ordinary'
    case_out.mkdir(parents=True, exist_ok=True)
    deck = case_out / case.deck.name
    prepare_deck(case.deck, deck)
    prefix = case_out / 'out'
    env = os.environ.copy()
    env.update({'MAD_NUM_THREADS': '2', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    command = [str(madqc), '--wf=response', str(deck), f'--prefix={prefix}']
    log = case_out / 'run.log'
    with log.open('w') as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=env, check=False)
    return {
        'run_id': f'{case.case_id}-ordinary', 'molecule': case.molecule, 'workflow': 'response',
        'source_branch': 'madqc/step4-astra', 'ground_state_endpoint': f'{case.case_id}-gs-{case.protocols[0]:.0e}',
        'protocol': {'thresholds': list(case.protocols)}, 'phase': 'response_iteration',
        'sampling_phase': 'ordinary', 'classification': 'ordinary',
        'artifacts': {'output_dir': str(case_out), 'deck': str(deck), 'log': str(log)},
        'status': 'completed' if completed.returncode == 0 else 'failed', 'returncode': completed.returncode,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--case', choices=[c.case_id for c in CASES], action='append')
    parser.add_argument('--madqc', type=Path, default=DEFAULT_MADQC)
    parser.add_argument('--output-root', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    selected = [c for c in CASES if not args.case or c.case_id in args.case]
    if not args.dry_run and not args.execute:
        parser.error('select --dry-run or --execute')
    args.output_root.mkdir(parents=True, exist_ok=True)
    if args.dry_run:
        for row in plan_records(args.output_root):
            print(json.dumps(row, sort_keys=True))
        return 0
    if not args.madqc.exists():
        parser.error(f'MADQC binary not found: {args.madqc}')
    results = [execute_case(case, args.madqc, args.output_root) for case in selected]
    manifest = args.output_root / 'run-results.json'
    manifest.write_text(json.dumps({'source': str(SOURCE), 'results': results}, indent=2) + '\n')
    print(manifest)
    for result in results:
        print(result['run_id'], result['status'], result['returncode'])
    return 0 if all(r['status'] == 'completed' for r in results) else 1

if __name__ == '__main__':
    raise SystemExit(main())
