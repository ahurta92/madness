"""Validation helpers for portable MADQC convergence evidence."""

from __future__ import annotations

import math
from collections.abc import Mapping
from numbers import Real

RUN_CLASSIFICATIONS = frozenset({"ordinary", "forced", "configured-progress"})


def _require(mapping: Mapping, key: str) -> object:
    if key not in mapping:
        raise ValueError(f"missing required field: {key}")
    return mapping[key]


def _finite_number(value: object, label: str) -> None:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(float(value)):
        raise ValueError(f"{label} must be a finite numeric value")


def validate_metric_record(record: Mapping) -> None:
    """Validate one scalar metric record.

    Unavailable values are represented as ``available=False`` and
    ``value=None``.  A numeric zero is a real measurement and must not be
    used as an unavailable sentinel.
    """
    if not isinstance(record, Mapping):
        raise ValueError("metric record must be a mapping")
    value = _require(record, "value")
    available = _require(record, "available")
    aggregation = _require(record, "aggregation")
    units = _require(record, "units")
    if not isinstance(available, bool):
        raise ValueError("metric available must be boolean")
    if not isinstance(aggregation, str) or not aggregation:
        raise ValueError("metric aggregation must be a non-empty string")
    if not isinstance(units, str) or not units:
        raise ValueError("metric units must be a non-empty string")
    if available:
        _finite_number(value, "available metric value")
    elif value is not None:
        raise ValueError("unavailable metric must use value=None")
    if "per_function" in record:
        values = record["per_function"]
        if values is not None:
            if not isinstance(values, list):
                raise ValueError("per_function must be a list or null")
            for index, item in enumerate(values):
                _finite_number(item, f"per_function[{index}]")


def validate_run_record(record: Mapping) -> None:
    """Validate a JSON-compatible run record."""
    if not isinstance(record, Mapping):
        raise ValueError("run record must be a mapping")
    required = (
        "run_id", "molecule", "workflow", "source_branch",
        "ground_state_endpoint", "protocol", "phase", "sampling_phase",
        "classification", "artifacts",
    )
    for key in required:
        _require(record, key)
    for key in ("run_id", "molecule", "workflow", "source_branch", "phase", "sampling_phase"):
        if not isinstance(record[key], str) or not record[key]:
            raise ValueError(f"{key} must be a non-empty string")
    classification = record["classification"]
    if classification not in RUN_CLASSIFICATIONS:
        raise ValueError(f"unknown classification: {classification}")
    if record["workflow"] == "response" and not record["ground_state_endpoint"]:
        raise ValueError("response records require ground_state_endpoint")
    if record["ground_state_endpoint"] is not None and not isinstance(record["ground_state_endpoint"], str):
        raise ValueError("ground_state_endpoint must be a string or null")
    if not isinstance(record["protocol"], Mapping):
        raise ValueError("protocol must be a mapping")
    if not isinstance(record["artifacts"], Mapping):
        raise ValueError("artifacts must be a mapping")
    metrics = record.get("metrics", {})
    if not isinstance(metrics, Mapping):
        raise ValueError("metrics must be a mapping")
    for name, metric in metrics.items():
        if not isinstance(name, str) or not name:
            raise ValueError("metric names must be non-empty strings")
        validate_metric_record(metric)
