#!/usr/bin/env python3
"""
Offline validator for property-component parallel scheduling.

Checks three artifacts:
1) claim locks:   property_component_claims/run_*.{beta,raman}.task*.lock
2) component shards: properties_components.group*.json
3) optional logs: response_console.group*.log (duplicate VBC writes)

This is intended for queue-down/offline debugging of the dynamic component
allocation path.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


CLAIM_FILE_RE = re.compile(
    r"^(?P<prefix>.+)\.(?P<kind>beta|raman)\.task(?P<task_id>\d+)\.lock$"
)
CLAIM_PAYLOAD_RE = re.compile(
    r"group=(?P<group>\S+)\s+freq_b=(?P<freq_b>[-+0-9.eE]+)\s+"
    r"freq_c=(?P<freq_c>[-+0-9.eE]+)\s+B=(?P<B>\S+)\s+C=(?P<C>\S+)"
)
SHARD_FILE_RE = re.compile(r"properties_components\.group(?P<gid>\d+)\.json$")
VBC_WRITE_RE = re.compile(r"Wrote VBC to\s+(?P<fname>\S+)")


def norm_freq(x: float, digits: int = 12) -> str:
    return f"{x:.{digits}g}"


@dataclass(frozen=True)
class ClaimTask:
    kind: str
    task_id: int
    group: str
    freq_b: float
    freq_c: float
    b_desc: str
    c_desc: str
    path: Path

    def key(self) -> Tuple[str, str, str, str]:
        # kind, freqB, freqC, B, C
        return (
            self.kind,
            norm_freq(self.freq_b),
            norm_freq(self.freq_c),
            self.b_desc,
            self.c_desc,
        )


@dataclass(frozen=True)
class ComponentRow:
    property_name: str
    component: Tuple[str, ...]
    freq_b: float
    freq_c: Optional[float]
    group_id: int
    file: Path

    def full_key(self) -> Tuple[str, Tuple[str, ...], str, str]:
        return (
            self.property_name,
            self.component,
            norm_freq(self.freq_b),
            norm_freq(self.freq_c if self.freq_c is not None else 0.0),
        )

    def claim_match_key(self) -> Optional[Tuple[str, str, str, str, str]]:
        if len(self.component) < 3:
            return None
        if self.property_name == "hyperpolarizability":
            kind = "beta"
        elif self.property_name == "raman":
            kind = "raman"
        else:
            return None
        if self.freq_c is None:
            return None
        return (
            kind,
            norm_freq(self.freq_b),
            norm_freq(self.freq_c),
            self.component[1],
            self.component[2],
        )


def find_claim_files(run_dir: Path, claim_glob: str) -> List[Path]:
    return sorted(p for p in run_dir.rglob(claim_glob) if p.is_file())


def parse_claim_task(path: Path) -> Optional[ClaimTask]:
    m = CLAIM_FILE_RE.match(path.name)
    if not m:
        return None
    payload = path.read_text(encoding="utf-8", errors="replace")
    pm = CLAIM_PAYLOAD_RE.search(payload)
    if not pm:
        return None
    return ClaimTask(
        kind=m.group("kind"),
        task_id=int(m.group("task_id")),
        group=pm.group("group"),
        freq_b=float(pm.group("freq_b")),
        freq_c=float(pm.group("freq_c")),
        b_desc=pm.group("B"),
        c_desc=pm.group("C"),
        path=path,
    )


def load_json_array(path: Path) -> List[dict]:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return []
    if not isinstance(obj, list):
        return []
    return [row for row in obj if isinstance(row, dict)]


def collect_component_rows(run_dir: Path, shard_glob: str) -> List[ComponentRow]:
    rows: List[ComponentRow] = []
    for path in sorted(run_dir.rglob(shard_glob)):
        m = SHARD_FILE_RE.match(path.name)
        if not m:
            continue
        gid = int(m.group("gid"))
        for row in load_json_array(path):
            prop = row.get("property")
            comp = row.get("component")
            freq_b = row.get("freqB")
            freq_c = row.get("freqC")
            if not isinstance(prop, str):
                continue
            if not isinstance(comp, list) or not all(isinstance(x, str) for x in comp):
                continue
            if not isinstance(freq_b, (int, float)):
                continue
            if freq_c is not None and not isinstance(freq_c, (int, float)):
                continue
            rows.append(
                ComponentRow(
                    property_name=prop,
                    component=tuple(comp),
                    freq_b=float(freq_b),
                    freq_c=float(freq_c) if freq_c is not None else None,
                    group_id=gid,
                    file=path,
                )
            )
    return rows


def collect_log_vbc_writes(run_dir: Path, log_glob: str) -> Dict[str, Set[Path]]:
    writes: Dict[str, Set[Path]] = defaultdict(set)
    for path in sorted(run_dir.rglob(log_glob)):
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        for line in text.splitlines():
            m = VBC_WRITE_RE.search(line)
            if m:
                writes[m.group("fname")].add(path)
    return writes


def summarize_claim_distribution(tasks: Iterable[ClaimTask]) -> Dict[str, Counter]:
    dist: Dict[str, Counter] = {"beta": Counter(), "raman": Counter()}
    for t in tasks:
        if t.kind in dist:
            dist[t.kind][t.group] += 1
    return dist


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Validate dynamic property-component scheduling artifacts."
    )
    ap.add_argument(
        "--run-dir",
        default=".",
        help="Root directory containing claim files/shards/logs (default: cwd).",
    )
    ap.add_argument(
        "--claim-glob",
        default="*.task*.lock",
        help="Glob for claim files (searched recursively).",
    )
    ap.add_argument(
        "--shard-glob",
        default="properties_components.group*.json",
        help="Glob for component shard files (searched recursively).",
    )
    ap.add_argument(
        "--log-glob",
        default="response_console.group*.log",
        help="Glob for subgroup logs (searched recursively).",
    )
    ap.add_argument(
        "--skip-log-check",
        action="store_true",
        help="Skip duplicate VBC-write scan in logs.",
    )
    ap.add_argument(
        "--tol",
        type=float,
        default=1e-10,
        help="Float tolerance for frequency comparisons.",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero if mismatches are found.",
    )
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.exists():
        print(f"run-dir does not exist: {run_dir}", file=sys.stderr)
        return 2

    claim_files = find_claim_files(run_dir, args.claim_glob)
    claim_tasks: List[ClaimTask] = []
    unparsable_claims: List[Path] = []
    for cf in claim_files:
        parsed = parse_claim_task(cf)
        if parsed is None:
            unparsable_claims.append(cf)
        else:
            claim_tasks.append(parsed)

    rows = collect_component_rows(run_dir, args.shard_glob)

    # Duplicate rows across groups (same property/component/frequencies).
    row_groups: Dict[Tuple[str, Tuple[str, ...], str, str], Set[int]] = defaultdict(set)
    for r in rows:
        row_groups[r.full_key()].add(r.group_id)
    duplicate_rows = {k: g for k, g in row_groups.items() if len(g) > 1}

    # Index rows by claim key (kind, freqB, freqC, B, C).
    rows_by_claim_key: Dict[Tuple[str, str, str, str, str], List[ComponentRow]] = (
        defaultdict(list)
    )
    for r in rows:
        mk = r.claim_match_key()
        if mk is not None:
            rows_by_claim_key[mk].append(r)

    # Expected A-components inferred from rows.
    expected_a: Dict[str, Set[str]] = {"beta": set(), "raman": set()}
    for r in rows:
        mk = r.claim_match_key()
        if mk is None:
            continue
        kind = mk[0]
        if r.component:
            expected_a[kind].add(r.component[0])

    missing_claim_outputs: List[ClaimTask] = []
    partial_claim_outputs: List[Tuple[ClaimTask, int, int]] = []
    for t in claim_tasks:
        key = (
            t.kind,
            norm_freq(t.freq_b),
            norm_freq(t.freq_c),
            t.b_desc,
            t.c_desc,
        )
        matched = rows_by_claim_key.get(key, [])
        if not matched:
            missing_claim_outputs.append(t)
            continue
        expected = len(expected_a[t.kind]) if expected_a[t.kind] else 0
        if expected > 0:
            got_a = {r.component[0] for r in matched if len(r.component) > 0}
            if len(got_a) < expected:
                partial_claim_outputs.append((t, len(got_a), expected))

    # Optional duplicate VBC writes from logs.
    duplicate_vbc_writes: Dict[str, Set[Path]] = {}
    if not args.skip_log_check:
        writes = collect_log_vbc_writes(run_dir, args.log_glob)
        duplicate_vbc_writes = {k: v for k, v in writes.items() if len(v) > 1}

    # Print summary.
    print("== Property Component Validator ==")
    print(f"run_dir: {run_dir}")
    print(f"claim_files: {len(claim_files)}")
    print(f"parsed_claim_tasks: {len(claim_tasks)}")
    print(f"component_rows: {len(rows)}")
    print(f"duplicate_row_keys_across_groups: {len(duplicate_rows)}")
    print(f"missing_claim_outputs: {len(missing_claim_outputs)}")
    print(f"partial_claim_outputs: {len(partial_claim_outputs)}")
    if not args.skip_log_check:
        print(f"duplicate_vbc_writes: {len(duplicate_vbc_writes)}")

    dist = summarize_claim_distribution(claim_tasks)
    for kind in ("beta", "raman"):
        if dist[kind]:
            parts = ", ".join(f"group {g}: {n}" for g, n in sorted(dist[kind].items()))
            print(f"{kind}_claim_distribution: {parts}")

    if unparsable_claims:
        print("\nUnparsable claim files:")
        for p in unparsable_claims[:20]:
            print(f"  - {p}")
        if len(unparsable_claims) > 20:
            print(f"  ... {len(unparsable_claims) - 20} more")

    if missing_claim_outputs:
        print("\nMissing claim outputs (no matching component rows):")
        for t in missing_claim_outputs[:20]:
            print(
                f"  - {t.kind} task={t.task_id} group={t.group} "
                f"freq_b={t.freq_b} freq_c={t.freq_c} B={t.b_desc} C={t.c_desc}"
            )
        if len(missing_claim_outputs) > 20:
            print(f"  ... {len(missing_claim_outputs) - 20} more")

    if partial_claim_outputs:
        print("\nPartial claim outputs (not all inferred A-components present):")
        for t, got, exp in partial_claim_outputs[:20]:
            print(
                f"  - {t.kind} task={t.task_id} group={t.group} "
                f"B={t.b_desc} C={t.c_desc} got_A={got}/{exp}"
            )
        if len(partial_claim_outputs) > 20:
            print(f"  ... {len(partial_claim_outputs) - 20} more")

    if duplicate_rows:
        print("\nDuplicate component rows across groups:")
        for key, gids in list(duplicate_rows.items())[:20]:
            prop, comp, fb, fc = key
            print(f"  - {prop} {comp} fb={fb} fc={fc} groups={sorted(gids)}")
        if len(duplicate_rows) > 20:
            print(f"  ... {len(duplicate_rows) - 20} more")

    if duplicate_vbc_writes:
        print("\nDuplicate VBC writes across logs:")
        for fname, files in list(duplicate_vbc_writes.items())[:20]:
            print(f"  - {fname}")
            for f in sorted(files):
                print(f"      {f}")
        if len(duplicate_vbc_writes) > 20:
            print(f"  ... {len(duplicate_vbc_writes) - 20} more")

    has_issues = (
        bool(unparsable_claims)
        or bool(missing_claim_outputs)
        or bool(partial_claim_outputs)
        or bool(duplicate_rows)
        or bool(duplicate_vbc_writes)
    )
    if args.strict and has_issues:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
