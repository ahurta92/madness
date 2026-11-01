"""Read a named run's records as the tests use them (D-105).

A stage's results.json keeps this run's judgement and names where the
generation's scientific report is (convergence_report.published_in); the
root results.json lists the stages by path. stage_record returns a stage's
record with convergence_report merged from the manifest it names, so a test
reads one dictionary as before; stages returns every stage's, in run order.
"""
import json
from pathlib import Path


def stage_record(stage_dir):
    stage_dir = Path(stage_dir)
    record = json.loads((stage_dir / 'results.json').read_text())
    report = record.get('convergence_report', {})
    published = report.get('published_in')
    if published:
        manifest = json.loads((stage_dir / published['manifest']).read_text())
        merged = manifest.get('metrics', {}).get('convergence_report', {})
        merged = dict(merged)
        merged.update(report)
        record['convergence_report'] = merged
    return record


def stages(root):
    root = Path(root)
    index = json.loads((root / 'results.json').read_text())
    return [stage_record((root / entry['path']).parent) for entry in index['stages']]
