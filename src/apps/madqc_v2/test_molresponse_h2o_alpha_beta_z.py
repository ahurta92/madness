#!/usr/bin/env python3

import argparse
import json
import os
import sys
import subprocess

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import cleanup, skip_on_small_machines


def response_rows_to_map(rows):
    mapped = {}
    for row in rows:
        key = (
            row["property"],
            tuple(row["component"]),
            round(float(row["freqB"]), 12),
            None if "freqC" not in row else round(float(row["freqC"]), 12),
        )
        mapped[key] = None if "value" not in row else float(row["value"])
    return mapped

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="H2O response regression test (alpha/beta, z direction)"
    )
    parser.add_argument(
        "--reference_directory",
        action="store",
        default="@SRCDIR@",
        help="directory containing test input and reference output",
    )
    args = parser.parse_args()

    # This response test is expensive on small machines.
    try:
        if skip_on_small_machines():
            print("Skipping this verylong test on small machines")
            sys.exit(77)
    except Exception:
        print("Unable to evaluate machine size from MAD_NUM_THREADS, skipping test")
        sys.exit(77)

    print("Testing @BINARY@/@TESTCASE@")
    print("reference files found in directory:", args.reference_directory)

    prefix = "mad_@BINARY@_@TESTCASE@"
    outputfile = prefix + ".calc_info.json"
    inputfile = os.path.join(args.reference_directory, "test_molresponse_h2o_alpha_beta_z.in")
    referencefile = os.path.join(args.reference_directory, prefix + ".calc_info.ref.json")

    if not os.path.exists(inputfile):
        print("Input file not found:", inputfile)
        sys.exit(1)
    if not os.path.exists(referencefile):
        print("Reference file not found:", referencefile)
        print("Skipping until expected output JSON is added.")
        sys.exit(77)

    cleanup(prefix)

    cmd = f"./@BINARY@ --wf=response --prefix={prefix} {inputfile}"
    print("executing\n", cmd)
    p = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    print("finished run")
    print(p.stdout)
    if p.stderr:
        print(p.stderr)
    exitcode = p.returncode
    print("exitcode", exitcode)
    if exitcode != 0:
        sys.exit(exitcode)

    with open(outputfile, "r", encoding="utf-8") as f:
        got = json.load(f)
    with open(referencefile, "r", encoding="utf-8") as f:
        ref = json.load(f)

    success = True

    # Basic structure checks
    if len(got["tasks"]) != len(ref["tasks"]):
        print("Mismatch in number of tasks:", len(got["tasks"]), len(ref["tasks"]))
        success = False

    if len(got["tasks"]) < 2 or len(ref["tasks"]) < 2:
        print("Expected at least two tasks (SCF + response) in output/reference.")
        sys.exit(1)

    # Ground-state checks
    got_scf = got["tasks"][0]
    ref_scf = ref["tasks"][0]
    if got_scf["model"] != ref_scf["model"]:
        print("SCF model mismatch:", got_scf["model"], ref_scf["model"])
        success = False
    if abs(float(got_scf["properties"]["energy"]) - float(ref_scf["properties"]["energy"])) > 1e-4:
        print(
            "SCF energy mismatch:",
            got_scf["properties"]["energy"],
            ref_scf["properties"]["energy"],
        )
        success = False

    # Response checks
    got_rsp = got["tasks"][1]
    ref_rsp = ref["tasks"][1]
    if got_rsp["type"] != ref_rsp["type"]:
        print("Response type mismatch:", got_rsp["type"], ref_rsp["type"])
        success = False

    got_rows = response_rows_to_map(got_rsp["properties"]["response_properties"])
    ref_rows = response_rows_to_map(ref_rsp["properties"]["response_properties"])

    if set(got_rows.keys()) != set(ref_rows.keys()):
        print("Response property key mismatch")
        print("Only in output:", sorted(set(got_rows.keys()) - set(ref_rows.keys())))
        print("Only in reference:", sorted(set(ref_rows.keys()) - set(got_rows.keys())))
        success = False
    else:
        for key in sorted(ref_rows.keys()):
            gval = got_rows[key]
            rval = ref_rows[key]
            if gval is None and rval is None:
                continue
            if gval is None or rval is None:
                print("Missing response value for key:", key, gval, rval)
                success = False
                continue
            if abs(gval - rval) > 1e-4:
                print("Response value mismatch for key:", key, gval, rval)
                success = False

    print("final success:", success)

    sys.exit(0 if success else 1)
