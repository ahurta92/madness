#!/usr/bin/env python3
"""Review regressions: plotting and recovery before a requested bundle exists."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--work")
    args = parser.parse_args()
    work = Path(args.work or tempfile.mkdtemp(prefix="scf-recovery-")).resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               OMPI_MCA_hwloc_base_binding_policy="none")
    env.setdefault("MAD_NUM_THREADS", "2")
    launch = [] if args.ranks == 1 else [args.mpiexec, "-n", str(args.ranks)]
    failures = []
    serial = 0

    def run(case, settings="", extra=()):
        nonlocal serial
        directory = work / case
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "input.in").write_text(f"""dft
 xc hf
 localize canon
 aobasis 6-31gss
 protocol [1e-4]
 convergence.criteria [bsh_residual]
 convergence.bsh_residual 1000
 maxiter 0
 print_level 0
 {settings}
end
molecule
 units atomic
 He 0 0 0
end
""")
        result = subprocess.run(launch + [str(Path(args.binary).resolve()), "input.in", "--prefix=run"] + list(extra),
                                cwd=directory, env=env, capture_output=True, text=True, timeout=180)
        serial += 1
        (work / f"{serial:02d}-{case}.log").write_text(result.stdout + result.stderr)
        return result.returncode

    def root(case):
        return work / case / "run"

    def block_main(case):
        task = root(case) / "checkpoints/ground_state/tasks/scf"
        task.mkdir(parents=True, exist_ok=True)
        if not (task / "latest").exists():
            (task / "latest").write_text('{"generation":null}')
        task.chmod(0o500)
        return task

    def preparation(case):
        marker = json.loads((root(case) / "checkpoints/scf_preparation.json").read_text())
        task = root(case) / "checkpoints" / marker["calculation_id"] / "tasks/scf"
        generation = json.loads((task / "latest").read_text())["generation"]
        return json.loads((task / generation / "manifest.json").read_text())

    assert run("plot") == 0
    if run("plot", "restart read_only\n plotcoul true\n plotcube true\n npt_plot 5") != 0:
        failures.append("read-only Coulomb plotting failed")

    settings = "nvalpha 1\n nv_extra 1\n nv_its 1"
    for case, extra, restart in [("continue", ["--continue"], "auto"), ("iterate", [], "iterate")]:
        task = block_main(case)
        try:
            assert run(case, settings) != 0
            assert preparation(case)["task"]["iteration"]["protocol_complete"]
        finally:
            task.chmod(0o700)
        if run(case, settings + f"\n restart {restart}", extra) != 0:
            failures.append(f"{case} rejected preparation-only progress")

    assert run("schedule") == 0
    task = block_main("schedule")
    try:
        assert run("schedule", settings) != 0
        original = preparation("schedule")
        assert run("schedule", settings + "\n nv_extra 2\n nv_step 1\n protocol [1e-6]") != 0
        changed = preparation("schedule")
        if changed["generation"] == original["generation"] or changed["task"]["definition"]["protocol"]["k"] != 8:
            failures.append("changed preparation schedule inherited obsolete completion")
    finally:
        task.chmod(0o700)
    print(json.dumps({"failures": failures, "logs": str(work)}))
    assert not failures, failures


if __name__ == "__main__":
    main()
