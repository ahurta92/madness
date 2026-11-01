#!/usr/bin/env python3
"""The from-store tier, level 1 (D-078, D-079): capture runs, pack them into
the fixture test_named_replay reads, and replay verdicts with no engine.

  replay_verdict.py capture <run-root> <deck> <capture-dir> [--label L]
                            [--since T] [--arg A ...]
      Capture a finished run: the deck, the run's id and arguments, and for
      every stage of <run-root> that this run wrote (its record names this
      run's attempt directory; the run's index results.json was written at
      or after T, the run's start in seconds since the epoch) its
      results.json and its latest published manifest, as
      <capture-dir>/<label>/<n>/ (n counts the captures of that label). A run
      that wrote no stage leaves nothing.

  replay_verdict.py pack <capture-dir> <fixture.json> [--commit C]
      Pack every captured stage into one fixture: the decks once each, and per
      stage the atom count the AppSpec fields take, the final manifest reduced
      to what the verdict reads (its state records and identity; a manifest
      without numerical.states, written before S4a, is kept whole), and the
      verdict the run recorded (results.json status completed = converged).

  replay_verdict.py run <capture-dir> --driver <test_named_replay>
      Pack to a temporary fixture and replay it.

A change that leaves the solver's inputs alone is checked against the fixture
in seconds; a change to what the solver reads still needs the solve tiers.
The fixture changes only with a decided record change.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# What shared_verdict reads of a state record.
STATE_KEYS = ("schema", "key", "kind", "identity", "protocol_step", "stop_reason",
              "diverged", "decided_by", "saved_measurements", "engine_measurements")


def capture(root, deck, out, label, since=0.0, args=None):
    root, out = Path(root), Path(out) / label
    index = root / "results.json"
    if not index.exists() or index.stat().st_mtime < since:
        return None  # this run wrote no record (refused at the deck check)
    run_id = json.loads(index.read_text())["run_id"]

    def written(stage):
        # The stages this run wrote: their record names this run's attempt.
        f = stage / "results.json"
        return (f.exists() and (stage / "checkpoints" / "latest").exists()
                and json.loads(f.read_text()).get("output_directory", "")
                .endswith("attempt-" + run_id))
    stages = [d for d in sorted(root.iterdir()) if d.is_dir()]
    if not any(written(d) for d in stages):
        return None
    out.mkdir(parents=True, exist_ok=True)
    n = len([d for d in out.iterdir() if d.is_dir()]) + 1
    dest = out / f"{n:03d}"
    dest.mkdir()
    shutil.copy(deck, dest / "deck.in")
    (dest / "run.json").write_text(json.dumps({"run_id": run_id, "args": args or []}))
    # Every stage's latest manifest (a response stage's atom count is read
    # from its reference's); results.json only for the stages this run wrote.
    for stage in stages:
        latest = stage / "checkpoints" / "latest"
        if not latest.exists():
            continue
        name = json.loads(latest.read_text())["generation"]
        s = dest / stage.name
        s.mkdir()
        shutil.copy(stage / "checkpoints" / name / "manifest.json", s / "manifest.json")
        if written(stage):
            shutil.copy(stage / "results.json", s / "results.json")
    return dest


def natom(case_dir, manifest):
    def symbols(m):
        return len(m["numerical"]["consumer_context"]["molecule"]["symbols"])
    if "consumer_context" in manifest.get("numerical", {}):
        return symbols(manifest)
    ref = (manifest.get("reference") or {}).get("stage")
    return symbols(json.loads((case_dir / ref / "manifest.json").read_text()))


def reduced(manifest):
    # Kept whole: a manifest written before S4a, whose records are rebuilt
    # from the rest, and one whose state list is empty (response's records
    # are then derived from its engine metadata).
    if not manifest["numerical"].get("states"):
        return manifest
    out = dict(manifest)
    out["numerical"] = {"states": [{k: s[k] for k in STATE_KEYS if k in s}
                                   for s in manifest["numerical"]["states"]]}
    out["position"], out["metrics"], out["targets"] = {}, {}, {}
    out["properties"], out["archives"] = [], []
    return out


def pack(captures, fixture, commit=None):
    decks, cases, seen = {}, [], set()
    for case_dir in sorted(Path(captures).glob("*/*")):
        # A run's records are captured once: a later capture naming the same
        # run is a refused run's view of its predecessor's records.
        run = json.loads((case_dir / "run.json").read_text())["run_id"]
        if run in seen:
            continue
        seen.add(run)
        deck = (case_dir / "deck.in").read_text()
        sha = hashlib.sha256(deck.encode()).hexdigest()[:16]
        decks[sha] = deck
        for stage in sorted(d for d in case_dir.iterdir() if d.is_dir()):
            if not (stage / "results.json").exists():
                continue
            record = json.loads((stage / "results.json").read_text())
            manifest = json.loads((stage / "manifest.json").read_text())
            cases.append({"id": f"{case_dir.parent.name}/{case_dir.name}/{stage.name}",
                          "deck": sha, "stage": stage.name,
                          "natom": natom(case_dir, manifest),
                          "manifest": reduced(manifest),
                          "recorded": record.get("status") == "completed"})
    data = {"schema": "madqc.replay.v1", "writer_commit": commit,
            "decks": decks, "cases": cases}
    Path(fixture).write_text(json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n")
    return len(cases), len(decks)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("capture")
    c.add_argument("root")
    c.add_argument("deck")
    c.add_argument("out")
    c.add_argument("--label", default="run")
    c.add_argument("--since", type=float, default=0.0)
    c.add_argument("--arg", action="append", default=[],
                   help="one of the run's arguments, recorded in run.json")
    p = sub.add_parser("pack")
    p.add_argument("captures")
    p.add_argument("fixture")
    p.add_argument("--commit")
    r = sub.add_parser("run")
    r.add_argument("captures")
    r.add_argument("--driver", required=True)
    a = ap.parse_args()
    if a.cmd == "capture":
        print(capture(a.root, a.deck, a.out, a.label, a.since, a.arg))
    elif a.cmd == "pack":
        print("%d cases, %d decks" % pack(a.captures, a.fixture, a.commit))
    else:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Path(tmp) / "fixture.json"
            pack(a.captures, fixture)
            sys.exit(subprocess.call([a.driver, str(fixture)]))


if __name__ == "__main__":
    main()
