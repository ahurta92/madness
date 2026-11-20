#!/usr/bin/env python3
"""Changed or newly iterated orbitals cannot inherit an old completion claim."""
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
    args = parser.parse_args()
    binary = str(Path(args.binary).resolve())
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               OMPI_MCA_hwloc_base_binding_policy="none")
    env.setdefault("MAD_NUM_THREADS", "2")
    launch = [] if args.ranks == 1 else [args.mpiexec, "-n", str(args.ranks)]
    with tempfile.TemporaryDirectory(prefix="workflow-requalification-") as directory:
        work = Path(directory)

        def run(engine, settings, *, prefix=None):
            prefix = prefix or engine
            convergence = ("econv 1e-4\n dconv 1e-3\n save true" if engine == "nemo" else
                           "convergence.criteria [bsh_residual, density_change, energy_change]")
            deck = f"""dft
 xc hf
 protocol [1e-4]
 k 6
 {convergence}
 {settings}
end
{'nemo' + chr(10) + 'end' if engine == 'nemo' else ''}
molecule
 units atomic
 He 0 0 0
end
"""
            (work / "input.in").write_text(deck)
            proc = subprocess.run(launch + [binary, f"--wf={engine}", "input.in",
                                           f"--prefix={prefix}"],
                                  cwd=work, env=env, capture_output=True, text=True,
                                  timeout=180)
            assert proc.returncode in (0, 1), proc.stdout + proc.stderr
            return proc

        # Nemo must discard loaded convergence stamps before changing functions.
        initial = run("nemo", "maxiter 30")
        assert initial.returncode == 0, initial.stdout + initial.stderr
        report = work / "nemo/task_0/nemo/nemo.calc_info.json"
        report.unlink()  # force the application through the native archive path
        stopped = run("nemo", "maxiter 1\n restart iterate")
        assert stopped.returncode == 1, "Nemo reused a pre-solve completion claim\n" + stopped.stdout
        record = json.loads(report.read_text())
        assert not record["outcome"]["complete"]
        assert record["convergence"]["converged_for_dconv"] > 1.0
        report.unlink()
        resumed = run("nemo", "maxiter 1\n restart auto")
        assert resumed.returncode == 1, "unfinished native archive was certified\n" + resumed.stdout

        # Reusing the summary cannot certify MOs after load_mos changes k.
        initial = run("scf", "maxiter 30\n prefix moldft", prefix="moldft")
        assert initial.returncode == 0, initial.stdout + initial.stderr
        report = work / "moldft/analysis/ground_state/moldft.calc_info.json"
        previous = report.read_bytes()
        readonly = run("scf", "maxiter 1\n prefix moldft\n xc lda\n restart read_only", prefix="moldft")
        assert readonly.returncode == 0, readonly.stdout + readonly.stderr
        assert report.read_bytes() == previous, "inspection must preserve the original report"
        # Inspection did not evaluate the changed request; the solver must do so.
        automatic = run("scf", "maxiter 1\n prefix moldft\n xc lda", prefix="moldft")
        assert automatic.returncode == 1, automatic.stdout + automatic.stderr
        restored = run("scf", "maxiter 30\n prefix moldft", prefix="moldft")
        assert restored.returncode == 0, restored.stdout + restored.stderr
        previous = report.read_bytes()
        # Shared native evidence must be checked after representation conversion.
        # One allowed update can now genuinely reconverge this near-solved He
        # state: SCF evaluates that update before stopping. With zero updates,
        # the missing comparison history must prevent an inherited engine pass.
        changed = run("scf", "maxiter 0\n prefix moldft\n k 7", prefix="moldft")
        assert changed.returncode == 1, "changed-k reload released unfinished SCF\n" + changed.stdout
        assert not json.loads(report.read_text())["outcome"]["complete"]

        # This prefix collides with the engine's auxiliary calc_info writer.
        # Even a real recompute must preserve the prior application report if
        # its checked temporary write fails.
        report.write_bytes(previous)
        (report.parent / "moldft.calc_info.json.tmp").mkdir()
        failed = run("scf", "maxiter 1\n prefix moldft\n k 8", prefix="moldft")
        assert failed.returncode == 1, failed.stdout + failed.stderr
        assert "Failed to write SCF checkpoint report" in failed.stdout
        assert report.read_bytes() == previous, "auxiliary writer overwrote the saved report"

        # A forced frozen-virtual solve must also invalidate the old archive's
        # completion before saving a capped result. Zero iterations is a
        # deterministic cap even when the loaded virtuals are already close.
        virtuals = run("scf", "maxiter 30\n nvalpha 0\n prefix virtuals", prefix="virtuals")
        assert virtuals.returncode == 0, virtuals.stdout + virtuals.stderr
        report = work / "virtuals/analysis/ground_state/moldft.calc_info.json"
        report.unlink()
        capped = run("scf", "maxiter 0\n nvalpha 1\n prefix virtuals\n freeze_occupied true\n restart iterate", prefix="virtuals")
        assert capped.returncode == 1, capped.stdout + capped.stderr
        record = json.loads(report.read_text())
        assert record["outcome"] == {"complete": False, "reason": "iteration_limit"}
        assert record["convergence"]["converged_for_thresh"] > 1.0
        assert record["convergence"]["converged_for_dconv"] > 1.0
        task = work / "virtuals/checkpoints/ground_state/tasks/virtuals"
        old_generation = json.loads((task / "latest").read_text())["generation"]
        # This coarse-grid case tests continuation, not the physical accuracy
        # of a helium virtual spectrum. Select an attainable explicit tolerance
        # and require actual updates plus measured evidence below that bound.
        resumed = run("scf", "maxiter 40\n nvalpha 1\n prefix virtuals\n freeze_occupied true\n convergence.bsh_residual 1e-2", prefix="virtuals")
        assert resumed.returncode == 0, resumed.stdout + resumed.stderr
        assert json.loads(report.read_text())["outcome"]["complete"]
        generation = json.loads((task / "latest").read_text())["generation"]
        assert generation != old_generation
        iteration = json.loads((task / generation / "manifest.json").read_text())["task"]["iteration"]
        assert 0 < iteration["updates_this_attempt"] <= 40
        residual = next(m for m in iteration["metrics"] if m["metric"] == "bsh_residual")
        assert residual["status"] == "measured" and residual["value"] <= 1e-2
    print("workflow requalification passed")


if __name__ == "__main__":
    main()
