import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from madqc_evidence_schema import validate_metric_record, validate_run_record


def base_run():
    return {
        "run_id": "h2-response-coarse-static",
        "molecule": "H2",
        "workflow": "response",
        "source_branch": "madqc/step4-astra",
        "ground_state_endpoint": "h2-gs-coarse",
        "protocol": {"threshold": 1e-4, "k": 6},
        "phase": "response_iteration",
        "sampling_phase": "ordinary",
        "classification": "ordinary",
        "artifacts": {"calc_info": "out.calc_info.json"},
    }


def test_response_requires_ground_state_endpoint():
    record = base_run()
    record["ground_state_endpoint"] = None
    with pytest.raises(ValueError, match="ground_state_endpoint"):
        validate_run_record(record)


def test_unknown_run_classification_is_rejected():
    record = base_run()
    record["classification"] = "calibrated_forced"
    with pytest.raises(ValueError, match="classification"):
        validate_run_record(record)


def test_unavailable_metric_cannot_have_numeric_zero():
    record = {
        "value": 0.0,
        "available": False,
        "aggregation": "max",
        "units": "hartree",
    }
    with pytest.raises(ValueError, match="unavailable"):
        validate_metric_record(record)


def test_open_shell_record_can_report_total_and_spin_density():
    record = base_run()
    record.update({"molecule": "LiH", "workflow": "moldft"})
    record["metrics"] = {
        "total_density_l2": {
            "value": 1.2e-6,
            "available": True,
            "aggregation": "max",
            "units": "l2",
        },
        "spin_density_l2": {
            "value": 3.4e-7,
            "available": True,
            "aggregation": "max",
            "units": "l2",
        },
    }
    assert validate_run_record(record) is None
