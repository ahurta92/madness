#!/usr/bin/env python3
"""The from-store tier, level 2 (D-079, D-080): a store of small published
generations that test_named_from_store.py continues from, so the read, the
plan, the publish, retention and the handoff run from a load, not a solve.

Never binaries in git: git holds test_store.json (this file's list of what
the store holds, with sizes and SHA-256) and this script. The files live at
<root>/<commit>/, where <root> is $MADQC_TEST_STORE (default
~/madqc-test-store) and <commit> the commit test_store.json names.

  make_test_store.py generate <madqc> [--mpiexec MPIEXEC]
      Run the store's decks into <root>/<commit>/ (<commit> = HEAD; refused
      when a C++ source or a deck under src/ differs from HEAD, or a deck is
      not in HEAD) and write test_store.json.
  make_test_store.py verify [<dir>]
      Every listed file present with its size and SHA-256; nothing else under
      a run's stage directories.
  make_test_store.py fetch
      Download $MADQC_TEST_STORE_URL/<commit>.tar (any URL urllib opens,
      file:// included), unpack it under <root> and verify.
  make_test_store.py archive <dir>
      Write <dir>/<commit>.tar from the store, for fetch.

The store regenerates only with a decided record change (the manifest or
state-record schema, the name form, what a generation holds); the commit that
regenerates it says so.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
LISTING = HERE / "test_store.json"
SOURCE = HERE.parent.parent.parent  # the repository root


def capped(deck):
    """The H2 ground state only, stopped at its iteration limit under written
    targets it cannot meet in three iterations (the reuse republish's input)."""
    gs = deck.split("molresponse")[0]
    return (gs.replace("moldft gs\n", "moldft gs\n  convergence.use_default_criteria false\n")
              .replace("  econv 1e-6\n", "  econv 1e-14\n")
              .replace("  maxiter 30\n", "  maxiter 3\n"))


# The store's runs: (name, deck, transform, ranks, the run is expected to end unconverged).
RUNS = [
    ("h2-r1", "named_h2_plumbing.in", None, 1, False),
    ("h2-r2", "named_h2_plumbing.in", None, 2, False),
    ("he-r1", "named_he_alpha.in", None, 1, False),
    ("he-r2", "named_he_alpha.in", None, 2, False),
    ("h2-capped-r1", "named_h2_plumbing.in", "capped", 1, True),
]
TRANSFORMS = {None: lambda d: d, "capped": capped}
RECORD_SCHEMA = {"manifest_version": 1, "state_schema": "madqc.state.v1",
                 "generation_names": "g-<stage>-<run>-<n>"}


def root():
    return Path(os.environ.get("MADQC_TEST_STORE", Path.home() / "madqc-test-store"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def listing():
    return json.loads(LISTING.read_text())


def store_dir(data=None):
    return root() / (data or listing())["commit"]


def deck_text(name, transform):
    return TRANSFORMS[transform]((HERE / name).read_text())


def generate(madqc, mpiexec):
    git = ["git", "-C", str(SOURCE)]
    changed = subprocess.run(git + ["diff", "--name-only", "HEAD", "--", "src"],
                             capture_output=True, text=True, check=True).stdout.split()
    # The C++ sources and the decks are the store's source: HEAD must hold
    # them as they are, so that the commit the listing names reproduces it.
    sources = [f for f in changed if f.endswith((".h", ".hpp", ".cc", ".cpp", ".in"))]
    sources += [f"src/apps/madqc/{deck}" for _, deck, _, _, _ in RUNS
                if subprocess.run(git + ["cat-file", "-e", f"HEAD:src/apps/madqc/{deck}"],
                                  capture_output=True).returncode != 0]
    if sources:
        sys.exit("refused: these sources differ from HEAD, so HEAD would not reproduce "
                 "the store: " + " ".join(sorted(set(sources))))
    commit = subprocess.run(git + ["rev-parse", "--short=9", "HEAD"],
                            capture_output=True, text=True, check=True).stdout.strip()
    out = root() / commit
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               OMPI_MCA_hwloc_base_binding_policy="none")
    env.setdefault("MAD_NUM_THREADS", "10")
    runs = []
    for name, deck, transform, ranks, unconverged in RUNS:
        text = deck_text(deck, transform)
        inp = out / f"{name}.in"
        inp.write_text(text)
        launch = [] if ranks == 1 else mpiexec.split() + ["-n", str(ranks)]
        r = subprocess.run(launch + [madqc, str(inp), f"--prefix={out / name}"], cwd=out,
                           env=env, capture_output=True, text=True)
        (out / f"{name}.log").write_text(r.stdout + r.stderr)
        if (r.returncode == 0) == unconverged:
            sys.exit(f"{name}: exit {r.returncode}; see {out / name}.log")
        files = [{"path": str(f.relative_to(out / name)), "bytes": f.stat().st_size,
                  "sha256": sha256(f)}
                 for f in sorted((out / name).rglob("*")) if f.is_file()]
        runs.append({"name": name, "deck": deck, "transform": transform,
                     "deck_sha256": hashlib.sha256(text.encode()).hexdigest(),
                     "ranks": ranks, "files": files})
        print(name, len(files), "files", sum(f["bytes"] for f in files), "bytes", flush=True)
    LISTING.write_text(json.dumps({"schema": "madqc.test-store.v1", "commit": commit,
                                   "record_schema": RECORD_SCHEMA, "runs": runs},
                                  indent=1) + "\n")
    print("wrote", out, "and", LISTING)


def verify(directory=None):
    data = listing()
    base = Path(directory) if directory else store_dir(data)
    problems = []
    for run in data["runs"]:
        listed = {f["path"]: f for f in run["files"]}
        top = base / run["name"]
        for path, f in listed.items():
            p = top / path
            if not p.is_file():
                problems.append(f"{run['name']}/{path}: missing")
            elif p.stat().st_size != f["bytes"] or sha256(p) != f["sha256"]:
                problems.append(f"{run['name']}/{path}: size or checksum differs")
        if top.is_dir():
            for p in top.rglob("*"):
                if p.is_file() and str(p.relative_to(top)) not in listed:
                    problems.append(f"{run['name']}/{p.relative_to(top)}: not in test_store.json")
    return problems


def fetch():
    data = listing()
    url = os.environ["MADQC_TEST_STORE_URL"].rstrip("/") + f"/{data['commit']}.tar"
    root().mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r, tarfile.open(fileobj=r, mode="r|") as tar:
        tar.extractall(root(), filter="data")
    return verify()


def archive(directory):
    data = listing()
    target = Path(directory) / f"{data['commit']}.tar"
    with tarfile.open(target, "w") as tar:
        for run in data["runs"]:
            for f in run["files"]:
                tar.add(store_dir(data) / run["name"] / f["path"],
                        arcname=f"{data['commit']}/{run['name']}/{f['path']}")
    return target


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("madqc")
    g.add_argument("--mpiexec", default="mpiexec --bind-to none")
    v = sub.add_parser("verify")
    v.add_argument("dir", nargs="?")
    sub.add_parser("fetch")
    a = sub.add_parser("archive")
    a.add_argument("dir")
    args = ap.parse_args()
    if args.cmd == "generate":
        generate(args.madqc, args.mpiexec)
    elif args.cmd in ("verify", "fetch"):
        problems = verify(args.dir) if args.cmd == "verify" else fetch()
        for p in problems[:40]:
            print(" ", p)
        print(f"{len(problems)} problems")
        sys.exit(1 if problems else 0)
    else:
        print(archive(args.dir))


if __name__ == "__main__":
    main()
