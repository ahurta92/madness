#!/usr/bin/env python3
"""Public SCF contract: reuse evidence, fresh checks, native integrity and inspection.

Numerical reference accuracy is covered separately by curated qctest decks.
This test deliberately uses loose/zero tolerances to distinguish controller
decisions without depending on the rate of scientific convergence.
"""
import argparse
import hashlib
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
    work = Path(args.work or tempfile.mkdtemp(prefix="shared-scf-")).resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               OMPI_MCA_hwloc_base_binding_policy="none")
    env.setdefault("MAD_NUM_THREADS", "2")
    launch = [] if args.ranks == 1 else [args.mpiexec, "-n", str(args.ranks)]
    counter = 0

    def run(directory, settings, *, success=True, extra=(), geometry="He 0.1 0 0"):
        nonlocal counter
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "input.in").write_text(f"""dft
 xc hf
 localize canon
 protocol [1e-4]
 convergence.bsh_residual 0
 maxiter 0
 print_level 0
 {settings}
end
molecule
 units atomic
 no_orient true
 {geometry}
end
""")
        command = launch + [str(Path(args.binary).resolve()), "input.in", "--wf=scf", "--prefix=run"]
        result = subprocess.run(command + list(extra), cwd=directory, env=env,
                                capture_output=True, text=True, timeout=240)
        counter += 1
        (work / f"{counter:02d}.log").write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == success, (command, result.returncode, result.stdout[-6000:], result.stderr)
        return directory / "run"

    def selected(root, task="scf"):
        directory = root / "checkpoints/ground_state/tasks" / task
        pointer = json.loads((directory / "latest").read_text())["generation"]
        # The shared store's pointer is a generation name, not scientific data.
        record = json.loads((directory / pointer / "manifest.json").read_text())
        return pointer, record

    def snapshot(root):
        return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in root.rglob("*") if p.is_file() and p.name != ".lock"}

    directory = work / "madqc"
    root = run(directory, "print_level 3", success=False)
    log = (work / "01.log").read_text()
    assert "Timing SCF" in log and "Preparation, evaluation and checks" in log
    assert "Checkpoint and reporting" in log and "Solver update" not in log
    assert "BSH residual, maximum" in log and "unmet" in log
    assert "Energy breakdown (hartree)" in log
    assert "Final orbital preparation: alpha" in log and "Transform orbitals" in log
    assert "Final orbital reevaluation" in log and "Final orbital spectrum: alpha" in log
    assert "Maximum spin-density change" in log
    assert "Kinetic change (Ha)" in log and "unavailable" in log
    assert "Checkpoint saved." in log and "No further update" in log
    assert "   conv " not in log, "legacy tolerance banner is misleading on the shared path"
    first, _ = selected(root)
    results = json.loads((root / "results.json").read_text())
    energy = results["calculations"][0]["properties"]["total_energy"]
    assert energy["status"] == "unavailable" and "value" not in energy
    (root / "analysis/ground_state/moldft.calc_info.json").write_text("broken analysis report")
    run(directory, "convergence.bsh_residual 1000\n print_level 1")
    log = (work / "02.log").read_text()
    assert "Timing SCF" not in log, "reuse must not report historical iteration work"
    assert "Wall time" in log
    assert "Checkpoint SCF" in log and "Load and checks" in log
    assert "Saved metrics against current tolerances" in log
    assert "no solver iterations required" in log
    assert "Checkpoint saved." not in log, "reuse does not publish"
    assert selected(root)[0] == first, "saved metrics should avoid another evaluation"
    run(directory, "maxiter 1\n print_level 1", success=False)
    log = (work / "03.log").read_text()
    assert log.count("Timing SCF") == 2
    assert log.count("Solver update") == 1
    assert "Continuing from saved functions with fresh convergence checks" in log
    assert "reported only" in log
    tightened, _ = selected(root)
    assert tightened != first, "a tighter request must iterate"
    before_model_change = snapshot(root / "checkpoints/ground_state/tasks/scf")
    run(directory, "convergence.bsh_residual 1000", geometry="He 0 0.1 0", success=False)
    assert snapshot(root / "checkpoints/ground_state/tasks/scf") == before_model_change
    moved = selected(root)[0]
    assert moved == tightened, "model rejection must preserve the saved checkpoint"
    run(directory, "protocol [1e-6]\n convergence.bsh_residual 1000",
        geometry="He 0.1 0 0")
    converted, record = selected(root)
    assert converted != moved
    definition = record["task"]["definition"]
    assert definition["protocol"] == {"threshold": 1e-6, "k": 8}
    before = snapshot(root / "checkpoints")
    run_record = (root / "checkpoints/run.json").read_bytes()
    results_record = (root / "results.json").read_bytes()
    run(directory, "restart read_only\n plotdens true\n plotcube true\n npt_plot 5")
    assert snapshot(root / "checkpoints") == before
    assert (root / "checkpoints/run.json").read_bytes() == run_record
    assert (root / "results.json").read_bytes() == results_record
    assert list(directory.glob("*.cube")), "inspection must load native functions for plotting"


    run(work / "absent", "", success=False, extra=["--continue"])
    assert not (work / "absent/run").exists(), "continue must not silently start a new calculation"
    run(work / "retired", "", success=False, extra=["--dft=dconv=1e-4", "--dft=maxiter=0"])
    assert not (work / "retired/run").exists()
    run(work / "invalid", "convergence.energy_change nan", success=False)
    assert not (work / "invalid/run").exists()
    json_directory = work / "json-invalid"
    json_directory.mkdir()
    (json_directory / "input.json").write_text(json.dumps({"dft": {"dconv": 1e-4}}))
    rejected = subprocess.run(launch + [str(Path(args.binary).resolve()), "input.json", "--dft=maxiter=0"],
                              cwd=json_directory, env=env, capture_output=True, text=True, timeout=60)
    (work / "json-invalid.log").write_text(rejected.stdout + rejected.stderr)
    assert rejected.returncode != 0 and "unknown input: dconv" in rejected.stdout + rejected.stderr
    # Reuse exact occupied evidence, extend the bundle and certify virtuals separately.
    directory = work / "madqc"
    root = run(directory, "protocol [1e-6]\n convergence.bsh_residual 1000\n freeze_occupied true\n nvalpha 1",
               geometry="He 0.1 0 0")
    virtual_generation, virtuals = selected(root, "virtuals")
    assert virtuals["task"]["inputs"][0]["generation"] == selected(root)[0]
    run(directory, "protocol [1e-6]\n freeze_occupied true\n nvalpha 1", success=False,
        geometry="He 0.1 0 0")
    assert selected(root, "virtuals")[0] == virtual_generation, "unqualified occupied source must block virtual iteration"
    prep = run(work / "preparation", "aobasis 6-31gss\n nvalpha 1\n nv_extra 1\n nv_its 1\n convergence.bsh_residual 1000")
    inventory = json.loads((prep / "checkpoints/run.json").read_text())
    assert any(c.get("retired") for c in inventory["calculations"]), "completed preparation must retire"
    assert len(list((prep / "checkpoints").rglob("scf.restartdata.00000"))) == 1, "preparation archives must be reclaimed"
    # Stop after preparation by denying publication in the main task only.
    directory = work / "interrupted-preparation"
    root = run(directory, "convergence.bsh_residual 1000")
    task_dir = root / "checkpoints/ground_state/tasks/scf"
    settings = "aobasis 6-31gss\n nvalpha 1\n nv_extra 1\n nv_its 1\n convergence.bsh_residual 1000"
    task_dir.chmod(0o500)
    try:
        run(directory, settings, success=False)
    finally:
        task_dir.chmod(0o700)
    progress = json.loads((root / "checkpoints/scf_preparation.json").read_text())
    assert progress.get("calculation_id"), "preparation progress must survive a failed main publication"
    run(directory, settings.replace("convergence.bsh_residual 1000", "convergence.bsh_residual 0") +
        "\n print_level 1", success=False)
    log = (work / f"{counter:02d}.log").read_text()
    preparation = log.split("SCF preparation restart:", 1)[1].split("Loaded saved preparation functions.", 1)[0]
    residual_row = next(line for line in preparation.splitlines() if "BSH residual, maximum" in line)
    assert "5.0000e-04" in residual_row, residual_row
    assert "Completed preparation bundle" in log
    timing = json.loads((root / "timings.json").read_text())["attempts"][-1]["timing"]
    preparation_loads = [node for node in timing["children"][0]["children"]
                         if node["name"] == "orbital_preparation" and
                         node.get("calculation_id") == progress["calculation_id"]]
    assert preparation_loads, "reused preparation needs its own timing identity"
    assert preparation_loads[0]["protocol"] == {"threshold": 1e-4, "k": 6}
    assert {node["name"] for node in preparation_loads[0]["children"]} == {"checkpoint_planning", "load"}
    assert all(node["seconds"] >= 0 for node in preparation_loads[0]["children"])
    run(directory, settings)
    assert json.loads((root / "checkpoints/scf_preparation.json").read_text()) == {"version": 1}
    assert len(list((root / "checkpoints").rglob("scf.restartdata.00000"))) == 1
    # Resuming a later intermediate protocol must show that protocol's floor,
    # even though workflow setup first visits the first requested protocol.
    directory = work / "later-protocol"
    run(directory, "", success=False)
    run(directory, "protocol [1e-3, 1e-4, 1e-6]\n print_level 1", success=False)
    log = (work / f"{counter:02d}.log").read_text()
    saved_table = log.split("Saved metrics against current tolerances", 1)[1].split(
        "Continuing from saved functions", 1)[0]
    residual_row = next(line for line in saved_table.splitlines() if "BSH residual, maximum" in line)
    assert "5.0000e-04" in residual_row, residual_row
    assert "5.0000e-03" not in residual_row
    # A selected native payload with changed bytes is an error, even when its metrics pass.
    root = work / "madqc/run"
    generation, _ = selected(root)
    archive = root / "checkpoints/ground_state/tasks/scf" / generation / "scf.restartdata.00000"
    with archive.open("ab") as stream:
        stream.write(b"corruption")
    run(work / "madqc", "convergence.bsh_residual 1000", success=False)
    print(f"shared SCF public contract passed; logs: {work}")


if __name__ == "__main__":
    main()
