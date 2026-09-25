#!/usr/bin/env python3
"""Assert the task-record envelope on a madqc calc_info.json.

Usage: check_task_envelope.py <prefix>.calc_info.json [--expect-response]
Exit 0 iff every assertion holds; prints one line per check.
"""
import json, sys

REQUIRED_TOP = ("schema_name", "schema_version", "provenance", "tasks")
REQUIRED_PROVENANCE = ("madness", "workflow", "hostname", "nproc", "threads")
REQUIRED_SCF = ("xc", "precision", "scf_iterations", "nuclear_repulsion_energy",
                "scf_one_electron_energy", "scf_two_electron_energy", "scf_kinetic_energy",
                "scf_nuclear_attraction_energy", "scf_coulomb_energy", "scf_total_energy")
REQUIRED_PRECISION = ("k", "thresh", "protocol", "econv", "dconv", "L", "ncoeff")

def main() -> int:
    path = sys.argv[1]
    expect_response = "--expect-response" in sys.argv[2:]
    ci = json.load(open(path))
    fails = 0
    def check(cond, label):
        nonlocal fails
        print(("  [PASS]  " if cond else "  [FAIL]  ") + label)
        if not cond: fails += 1
    for k in REQUIRED_TOP: check(k in ci, f"top-level '{k}' present")
    for k in REQUIRED_PROVENANCE: check(k in ci.get("provenance", {}), f"provenance.{k} present")
    check("git_commit" in ci.get("provenance", {}).get("madness", {}), "provenance.madness.git_commit present")
    scf = next((t for t in ci["tasks"] if t.get("type") in ("scf", "nemo")), None)
    check(scf is not None, "an SCF task entry has type scf|nemo")
    if scf is not None:
        for k in REQUIRED_SCF: check(k in scf["scf"], f"scf.{k} present")
        for k in REQUIRED_PRECISION: check(k in scf["scf"].get("precision", {}), f"scf.precision.{k} present")
        check(scf.get("precision") == scf["scf"].get("precision"), "precision mirrored at task top level")
        check(scf["convergence"].get("status") in ("converged", "unconverged"), "convergence.status set")
        check(scf["convergence"].get("iterations", -1) >= 1, "convergence.iterations >= 1")
        check("wall_s" in scf.get("provenance", {}), "task provenance.wall_s present")
        e = scf["scf"]
        # one- + two-electron + nuclear repulsion (+ pcm, disp) reproduces the total to the recorded precision
        need = ("scf_one_electron_energy", "scf_two_electron_energy", "nuclear_repulsion_energy", "scf_total_energy")
        if all(k in e for k in need):
            total = (e["scf_one_electron_energy"] + e["scf_two_electron_energy"] + e["nuclear_repulsion_energy"]
                     + e.get("scf_pcm_energy", 0.0) + e.get("scf_dispersion_correction_energy", 0.0))
            check(abs(total - e["scf_total_energy"]) < 1e-8, f"energy components sum to scf_total_energy (diff {total - e['scf_total_energy']:.2e})")
        else:
            check(False, "energy components present for the sum check")
    if expect_response:
        resp = next((t for t in ci["tasks"] if t.get("type") == "response"), None)
        check(resp is not None, "a response task entry exists")
        if resp is not None:
            check("wall_s" in resp.get("provenance", {}), "response task provenance.wall_s present")
    print(f"{fails} failure(s)")
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main())
