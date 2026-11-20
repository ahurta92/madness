#!/usr/bin/env python3
"""An unfinished SCF must not release response; saved work can finish on restart."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpiexec")
    args = parser.parse_args()
    binary = str(Path(args.binary).resolve())
    deck = (Path(args.source) / "src/examples/qc/response_he_alpha/response_he_alpha.in").read_text()
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               OMPI_MCA_hwloc_base_binding_policy="none")
    env.setdefault("MAD_NUM_THREADS", "2")
    launch = [] if args.ranks == 1 else [args.mpiexec, "-n", str(args.ranks)]
    with tempfile.TemporaryDirectory(prefix="workflow-completion-") as directory:
        work = Path(directory)

        def run(text):
            (work / "input.in").write_text(text)
            proc = subprocess.run(
                launch + [binary, "--wf=response", "input.in", "--prefix=run"],
                cwd=work, env=env, capture_output=True, text=True, timeout=180,
            )
            if proc.returncode not in (0, 1):
                raise AssertionError(proc.stdout + proc.stderr)
            return proc

        stopped = run(deck.replace("maxiter   20", "maxiter   1"))
        assert stopped.returncode == 1, "an unfinished workflow must return nonzero"
        assert not (work / "run/analysis/response").exists(), "response ran on an unfinished SCF"
        records = json.loads((work / "run.calc_info.json").read_text())["tasks"]
        assert len(records) == 1 and records[0].get("type") != "task_failed"
        assert records[0]["outcome"] == {"complete": False, "reason": "iteration_limit"}
        state_dir = work / "run"
        assert (state_dir / "analysis/ground_state/moldft.calc_info.json").is_file()
        assert list((state_dir / "checkpoints").rglob("*.restartdata.00000")), "unfinished orbitals were not saved"
        assert "caught an exception" not in stopped.stdout

        completed = run(deck)
        assert completed.returncode == 0, completed.stdout + completed.stderr
        records = json.loads((work / "run.calc_info.json").read_text())["tasks"]
        assert len(records) == 2 and records[0]["outcome"]["complete"]
        assert records[1]["type"] == "response"
        assert (work / "run/analysis/response").is_dir()

        # A report-write failure is an error, not a controlled solver stop;
        # all ranks must exit and the previously published SCF report survives.
        report = state_dir / "analysis/ground_state/moldft.calc_info.json"
        previous_report = report.read_bytes()
        (state_dir / "analysis/ground_state/moldft.calc_info.json.tmp").mkdir()
        failed = run(deck)
        assert failed.returncode == 1, failed.stdout + failed.stderr
        assert "Failed to write SCF checkpoint report" in failed.stdout
        assert report.read_bytes() == previous_report
        records = json.loads((work / "run.calc_info.json").read_text())["tasks"]
        assert len(records) == 1 and records[0]["type"] == "task_failed"
    print("workflow completion and restart passed")


if __name__ == "__main__":
    main()
