#!/usr/bin/env python3
"""Unit tests for run_qctest.compare: the tol / rtol / max / allow_zero semantics
every qctest check.json relies on (and min), the key-path selectors, "count", and the
per-check "file".

Run by hand:   python3 bin/test_run_qctest.py
Via ctest:     madness/test/qc/run_qctest_selftest/run   (labels qctest;short)
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_qctest  # noqa: E402


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.out, self.ref = d / "out.json", d / "ref.json"

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, out, ref):
        self.out.write_text(json.dumps(out))
        self.ref.write_text(json.dumps(ref))

    def compare(self, checks):
        return run_qctest.compare(self.out, self.ref, checks)

    def test_tol_passes_within_and_fails_outside(self):
        self.write({"e": 1.0005}, {"e": 1.0})
        self.assertTrue(self.compare([{"key": ["e"], "tol": 1e-3}]))
        self.assertFalse(self.compare([{"key": ["e"], "tol": 1e-4}]))

    def test_rtol_is_relative_to_the_reference(self):
        self.write({"b": 101.0}, {"b": 100.0})
        self.assertTrue(self.compare([{"key": ["b"], "rtol": 0.02}]))
        self.assertFalse(self.compare([{"key": ["b"], "rtol": 0.005}]))

    def test_rtol_rejects_a_zero_reference_unless_allow_zero(self):
        self.write({"b": 0.0}, {"b": 0.0})
        self.assertFalse(self.compare([{"key": ["b"], "rtol": 0.01}]))
        self.assertTrue(self.compare([{"key": ["b"], "rtol": 0.01, "allow_zero": True}]))

    def test_max_bounds_the_output_only(self):
        self.write({"it": 9}, {"it": 7})  # an exact int tol would fail; max must not
        self.assertTrue(self.compare([{"key": ["it"], "max": 12}]))
        self.assertFalse(self.compare([{"key": ["it"], "max": 8}]))

    def test_min_bounds_the_output_only(self):
        self.write({"ov": 0.998}, {"ov": 0.9995})  # an exact tol would fail; min must not
        self.assertTrue(self.compare([{"key": ["ov"], "min": 0.99}]))
        self.assertTrue(self.compare([{"key": ["ov"], "min": 0.998}]))  # inclusive
        self.assertFalse(self.compare([{"key": ["ov"], "min": 0.999}]))

    def test_min_does_not_read_the_reference(self):
        self.write({"ov": 0.998}, {})
        self.assertTrue(self.compare([{"key": ["ov"], "min": 0.99}]))
        self.ref.unlink()
        self.assertTrue(self.compare([{"key": ["ov"], "min": 0.99}]))

    def test_min_rejects_a_non_numeric_or_missing_value(self):
        self.write({"s": "x"}, {"s": "x"})
        self.assertFalse(self.compare([{"key": ["s"], "min": 0}]))
        self.assertFalse(self.compare([{"key": ["zz"], "min": 0}]))

    def test_min_and_max_together_bracket_the_value(self):
        self.write({"v": 5}, {"v": 5})
        self.assertTrue(self.compare([{"key": ["v"], "min": 4, "max": 6}]))
        self.assertFalse(self.compare([{"key": ["v"], "min": 6, "max": 9}]))
        self.assertFalse(self.compare([{"key": ["v"], "min": 1, "max": 4}]))

    def test_min_and_tol_both_apply(self):
        self.write({"v": 1.0}, {"v": 1.5})
        self.assertTrue(self.compare([{"key": ["v"], "min": 0.5, "tol": 0.6}]))
        self.assertFalse(self.compare([{"key": ["v"], "min": 0.5, "tol": 0.1}]))
        self.assertFalse(self.compare([{"key": ["v"], "min": 1.2, "tol": 0.6}]))

    def test_an_unknown_field_fails_the_check(self):
        self.write({"v": 1.0}, {"v": 1.0})
        self.assertFalse(self.compare([{"key": ["v"], "mni": 0.5}]))
        self.assertFalse(self.compare([{"key": ["v"], "tol": 0.1, "atol": 0.1}]))

    def test_max_rejects_a_non_numeric_value(self):
        self.write({"s": "converged"}, {"s": "converged"})
        self.assertFalse(self.compare([{"key": ["s"], "max": 1}]))

    def test_max_and_tol_both_apply(self):
        self.write({"it": 9}, {"it": 9})
        self.assertTrue(self.compare([{"key": ["it"], "max": 12, "tol": 0}]))
        self.write({"it": 9}, {"it": 8})
        self.assertFalse(self.compare([{"key": ["it"], "max": 12, "tol": 0}]))

    def test_missing_key_fails_for_every_kind(self):
        self.write({"a": 1.0}, {"a": 1.0})
        self.assertFalse(self.compare([{"key": ["zz"], "max": 1}]))
        self.assertFalse(self.compare([{"key": ["zz"], "rtol": 0.1}]))
        self.assertFalse(self.compare([{"key": ["zz"], "tol": 0.1}]))

    def test_strings_and_bools_still_compare_exactly(self):
        self.write({"s": "converged", "c": True}, {"s": "converged", "c": True})
        self.assertTrue(self.compare([{"key": ["s"], "tol": 0}, {"key": ["c"], "tol": 0}]))
        self.write({"s": "unconverged", "c": False}, {"s": "converged", "c": True})
        self.assertFalse(self.compare([{"key": ["s"], "tol": 0}]))
        self.assertFalse(self.compare([{"key": ["c"], "tol": 0}]))

    def test_zero_reference_with_tol_needs_allow_zero(self):
        self.write({"w": 0.0}, {"w": 0.0})
        self.assertFalse(self.compare([{"key": ["w"], "tol": 0}]))
        self.assertTrue(self.compare([{"key": ["w"], "tol": 0, "allow_zero": True}]))


class SelectorTests(unittest.TestCase):
    """Key paths that select a list item by a field's value, not by position."""

    RESULTS = {"calculations": [
        {"calculation_id": "ground_state", "properties": {"total_energy": {"value": -1.13}}},
        {"calculation_id": "response", "properties": {
            "polarizability": {"values": [
                {"frequency": 0.0, "status": "complete", "components": {"zz": 6.43}},
                {"frequency": 0.05, "status": "complete", "components": {"zz": 6.51}}]},
            "hyperpolarizability": {"values": [
                {"frequency_b": 0.0, "frequency_c": 0.0, "components": {"zzz": 315.7}},
                {"frequency_b": 0.0, "frequency_c": 0.05, "components": {"zzz": 320.1}}]},
            "excited_states": {"values": [
                {"root": "es_root_0000", "excitation_energy": 0.4658},
                {"root": "es_root_0001", "excitation_energy": 0.5}]},
            "two": {"values": [{"root": "r", "x": 1}, {"root": "r", "x": 2}]}}}]}

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.out, self.ref = d / "out.json", d / "ref.json"

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, out, ref):
        self.out.write_text(json.dumps(out))
        self.ref.write_text(json.dumps(ref))

    def response(self, *path):
        return ["calculations", {"calculation_id": "response"}, "properties", *path]

    def test_lookup_selects_by_string_field(self):
        self.assertEqual(run_qctest.lookup(self.RESULTS, self.response(
            "excited_states", "values", {"root": "es_root_0001"}, "excitation_energy")), 0.5)

    def test_lookup_selects_by_number_within_1e_9(self):
        path = self.response("polarizability", "values", {"frequency": 0.05 + 5e-10}, "components", "zz")
        self.assertEqual(run_qctest.lookup(self.RESULTS, path), 6.51)

    def test_lookup_number_outside_1e_9_matches_nothing(self):
        path = self.response("polarizability", "values", {"frequency": 0.05 + 1e-8}, "components", "zz")
        with self.assertRaisesRegex(run_qctest.SelectorError, "matched 0 of 2"):
            run_qctest.lookup(self.RESULTS, path)

    def test_lookup_selects_by_a_pair_of_fields(self):
        path = self.response("hyperpolarizability", "values",
                             {"frequency_b": 0.0, "frequency_c": 0.05}, "components", "zzz")
        self.assertEqual(run_qctest.lookup(self.RESULTS, path), 320.1)

    def test_lookup_two_matches_is_an_error(self):
        with self.assertRaisesRegex(run_qctest.SelectorError, "matched 2 of 2"):
            run_qctest.lookup(self.RESULTS, self.response("two", "values", {"root": "r"}))

    def test_lookup_selector_on_an_object_is_an_error(self):
        with self.assertRaises(TypeError):
            run_qctest.lookup(self.RESULTS, ["calculations", 0, {"calculation_id": "x"}])

    def test_selector_does_not_match_a_bool_to_a_number(self):
        self.assertFalse(run_qctest._matches({"converged": True}, {"converged": 1}))

    def test_compare_passes_with_selectors(self):
        self.write(self.RESULTS, self.RESULTS)
        self.assertTrue(self.compare([
            {"key": self.response("polarizability", "values", {"frequency": 0.0}, "components", "zz"),
             "tol": 1e-3},
            {"key": self.response("polarizability", "values", {"frequency": 0.0}, "status"), "tol": 0},
            {"key": ["calculations", {"calculation_id": "ground_state"}, "properties",
                     "total_energy", "value"], "tol": 1e-5}]))

    def compare(self, checks):
        return run_qctest.compare(self.out, self.ref, checks)

    def test_compare_fails_a_selector_with_no_match(self):
        self.write(self.RESULTS, self.RESULTS)
        key = self.response("excited_states", "values", {"root": "es_root_0009"}, "excitation_energy")
        self.assertFalse(self.compare([{"key": key, "tol": 1e-3}]))
        self.assertFalse(self.compare([{"key": key, "rtol": 1e-3}]))
        self.assertFalse(self.compare([{"key": key, "max": 1}]))

    def test_compare_fails_a_selector_with_two_matches(self):
        self.write(self.RESULTS, self.RESULTS)
        key = self.response("two", "values", {"root": "r"}, "x")
        self.assertFalse(self.compare([{"key": key, "tol": 0}]))
        self.assertFalse(self.compare([{"key": key, "max": 5}]))

    def test_compare_fails_when_only_the_output_lacks_the_item(self):
        out = json.loads(json.dumps(self.RESULTS))
        out["calculations"][1]["properties"]["polarizability"]["values"].pop()
        self.write(out, self.RESULTS)
        key = self.response("polarizability", "values", {"frequency": 0.05}, "components", "zz")
        self.assertFalse(self.compare([{"key": key, "tol": 1e-3}]))

    def test_selected_value_outside_tolerance_fails(self):
        out = json.loads(json.dumps(self.RESULTS))
        out["calculations"][1]["properties"]["polarizability"]["values"][1]["components"]["zz"] = 6.6
        self.write(out, self.RESULTS)
        key = self.response("polarizability", "values", {"frequency": 0.05}, "components", "zz")
        self.assertFalse(self.compare([{"key": key, "tol": 1e-3}]))
        self.assertTrue(self.compare([{"key": key, "tol": 0.1}]))

    def test_count_compares_list_lengths(self):
        key = self.response("excited_states", "values")
        ref = json.loads(json.dumps(self.RESULTS))
        self.write(self.RESULTS, ref)
        self.assertTrue(self.compare([{"key": key, "count": True, "tol": 0}]))
        ref["calculations"][1]["properties"]["excited_states"]["values"].pop()
        self.write(self.RESULTS, ref)
        self.assertFalse(self.compare([{"key": key, "count": True, "tol": 0}]))

    def test_count_with_max_and_on_a_non_list(self):
        self.write(self.RESULTS, self.RESULTS)
        key = self.response("excited_states", "values")
        self.assertTrue(self.compare([{"key": key, "count": True, "max": 2}]))
        self.assertFalse(self.compare([{"key": key, "count": True, "max": 1}]))
        self.assertFalse(self.compare([{"key": key + [0], "count": True, "tol": 0}]))


class FileSelectionTests(unittest.TestCase):
    """Each check reads the file its "file" names; calc_info.json by default."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.paths = {kind: (d / f"out.{kind}.json", d / f"ref.{kind}.json")
                      for kind in run_qctest.FILES}
        # The same key holds different values in the two files, so a check that
        # reads the wrong file fails.
        self.write("calc_info", {"e": 1.0}, {"e": 1.0})
        self.write("results", {"e": 2.0}, {"e": 2.0})

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, kind, out, ref):
        self.paths[kind][0].write_text(json.dumps(out))
        self.paths[kind][1].write_text(json.dumps(ref))

    def compare(self, checks, kinds=run_qctest.FILES):
        return run_qctest.compare({k: self.paths[k][0] for k in kinds},
                                  {k: self.paths[k][1] for k in kinds}, checks)

    def test_default_reads_calc_info(self):
        self.write("results", {"e": 2.0}, {"e": 9.0})
        self.assertTrue(self.compare([{"key": ["e"], "tol": 0}]))
        self.assertTrue(self.compare([{"key": ["e"], "tol": 0, "file": "calc_info"}]))

    def test_results_check_reads_results(self):
        self.write("calc_info", {"e": 1.0}, {"e": 9.0})
        self.assertTrue(self.compare([{"key": ["e"], "tol": 0, "file": "results"}]))
        self.write("results", {"e": 2.0}, {"e": 2.5})
        self.assertFalse(self.compare([{"key": ["e"], "tol": 0, "file": "results"}]))

    def test_mixed_case_reads_both(self):
        self.assertTrue(self.compare([{"key": ["e"], "tol": 0},
                                      {"key": ["e"], "max": 2.0, "file": "results"}]))
        self.assertFalse(self.compare([{"key": ["e"], "tol": 0},
                                       {"key": ["e"], "max": 1.5, "file": "results"}]))

    def test_unknown_file_fails(self):
        self.assertFalse(self.compare([{"key": ["e"], "tol": 0, "file": "metadata"}]))

    def test_a_file_not_given_fails_its_checks(self):
        self.assertFalse(self.compare([{"key": ["e"], "tol": 0, "file": "results"}],
                                      kinds=("calc_info",)))

    def test_a_missing_file_fails_its_checks(self):
        self.paths["results"][0].unlink()
        self.assertFalse(self.compare([{"key": ["e"], "tol": 0, "file": "results"}]))
        self.assertTrue(self.compare([{"key": ["e"], "tol": 0}]))

    def test_max_only_does_not_need_the_reference(self):
        self.paths["results"][1].unlink()
        self.assertTrue(self.compare([{"key": ["e"], "max": 3, "file": "results"}]))

    def test_file_kinds_and_paths(self):
        self.assertEqual(run_qctest.file_kinds([{"key": []}]), ["calc_info"])
        self.assertEqual(run_qctest.file_kinds([{"key": [], "file": "results"}]), ["results"])
        self.assertEqual(run_qctest.file_kinds([{"key": [], "file": "results"}, {"key": []}]),
                         ["calc_info", "results"])
        case, work = Path("/src/qc/case_x"), Path("/scratch")
        files = run_qctest.output_files(case, work, {})
        self.assertEqual(files["calc_info"], (work / "case_x.calc_info.json",
                                              case / "reference" / "case_x.calc_info.json"))
        self.assertEqual(files["results"], (work / "case_x" / "results.json",
                                            case / "reference" / "case_x.results.json"))

    def test_update_writes_one_reference_per_file(self):
        case = Path(self.tmp.name) / "case"
        case.mkdir()
        report = Path(self.tmp.name) / "case.out"
        report.write_text("log")
        run_qctest.update_reference(case, [self.paths["calc_info"][0], self.paths["results"][0]],
                                    report, ["case.calc_info.json", "case.results.json"])
        self.assertEqual(json.loads((case / "reference" / "case.results.json").read_text()), {"e": 2.0})
        self.assertEqual(json.loads((case / "reference" / "case.calc_info.json").read_text()), {"e": 1.0})
        self.assertTrue((case / "reference" / "case.out").is_file())


class ExpectErrorTests(unittest.TestCase):
    """A refusal case: the deck must be rejected, loudly, for the stated reason."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.out = self.dir / "case.calc_info.json"

    def tearDown(self):
        self.tmp.cleanup()

    def write_failed(self, error):
        self.out.write_text(json.dumps({"tasks": [
            {"type": "scf"},
            {"type": "task_failed", "task_index": 1, "error": error}]}))

    def check(self, exitcode, expected):
        return run_qctest.check_expected_error(self.out, exitcode, expected)

    def test_refusal_with_the_expected_reason_passes(self):
        self.write_failed("response: beta.or is not implemented (no signed-frequency source)")
        self.assertTrue(self.check(1, "beta.or is not implemented"))

    def test_a_run_that_succeeds_fails_the_case(self):
        self.out.write_text(json.dumps({"tasks": [{"type": "scf"}, {"type": "response"}]}))
        self.assertFalse(self.check(0, "beta.or is not implemented"))

    def test_a_failure_for_another_reason_fails_the_case(self):
        self.write_failed("eval: coordinate upper-bound error in dimension")
        self.assertFalse(self.check(1, "beta.or is not implemented"))

    def test_a_nonzero_exit_without_a_record_fails_the_case(self):
        self.assertFalse(self.check(1, "beta.or is not implemented"))

    def test_load_check_accepts_expect_error_without_checks(self):
        (self.dir / "check.json").write_text(json.dumps({"expect_error": "x"}))
        self.assertEqual(run_qctest.load_check(self.dir, False)["expect_error"], "x")


class ExpectFilesTests(unittest.TestCase):
    """Side outputs a case must produce (e.g. the ES analysis JSON)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_every_pattern_must_match_a_file(self):
        (self.dir / "a" / "b").mkdir(parents=True)
        (self.dir / "a" / "b" / "es_analysis__1e-04_k6.json").write_text("{}")
        self.assertTrue(run_qctest.check_expected_files(
            self.dir, ["a/*/es_analysis__*.json"]))
        self.assertFalse(run_qctest.check_expected_files(
            self.dir, ["a/*/es_analysis__*.json", "a/*/missing_*.json"]))

    def test_no_patterns_passes(self):
        self.assertTrue(run_qctest.check_expected_files(self.dir, []))


class ExternalReferenceTests(unittest.TestCase):
    """A reference holding another code's numbers (reference_source) must not be
    overwritten by --update with MADNESS's own result."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.case = Path(self.tmp.name) / "case"
        (self.case / "reference").mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def test_declared_source_is_reported(self):
        (self.case / "reference" / "case.calc_info.json").write_text(
            json.dumps({"reference_source": "DALTON t-aug-cc-pVQZ", "tasks": []}))
        self.assertEqual(run_qctest.external_reference_source(self.case),
                         "DALTON t-aug-cc-pVQZ")

    def test_ordinary_reference_is_not_external(self):
        (self.case / "reference" / "case.calc_info.json").write_text(
            json.dumps({"tasks": []}))
        self.assertIsNone(run_qctest.external_reference_source(self.case))

    def test_declared_source_in_a_results_reference_is_reported(self):
        (self.case / "reference" / "case.results.json").write_text(
            json.dumps({"reference_source": "DALTON t-aug-cc-pVQZ", "calculations": []}))
        self.assertEqual(run_qctest.external_reference_source(self.case),
                         "DALTON t-aug-cc-pVQZ")

    def test_no_reference_yet_is_not_external(self):
        self.assertIsNone(run_qctest.external_reference_source(self.case))


if __name__ == "__main__":
    unittest.main()
