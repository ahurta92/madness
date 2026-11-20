"""Shared response result checks independent of iteration counts and JSON formatting."""
from copy import deepcopy
from math import isclose


def assert_same_alpha_properties(before, after):
    """Native reuse preserves sources exactly; contractions may differ by roundoff."""
    expected = deepcopy(before)
    old_rows = expected['polarizability']['values']
    new_rows = after['polarizability']['values']
    assert len(old_rows) == len(new_rows)
    for old, new in zip(old_rows, new_rows):
        assert old['components'].keys() == new['components'].keys()
        for component, value in old['components'].items():
            assert isclose(value, new['components'][component], rel_tol=1e-12, abs_tol=1e-11)
        old['components'] = new['components']
    # Frequencies, units, availability, status and all other fields stay exact.
    assert expected == after


def assert_response_setup_timed(attempt):
    """The application interval includes reference/input preparation before tasks."""
    run = attempt['timing']
    response = next(n for n in run['children'] if n['name'] == 'response')
    setup = next(n for n in response['children'] if n['name'] == 'reference_and_input_setup')
    assert 0 <= setup['seconds'] <= response['seconds'] <= run['seconds']


def assert_energy_only_es_properties(properties):
    """An explicitly truncated diagnostic must not publish uncomputed zeros."""
    assert properties['excitation_energy']['components']
    assert not properties['oscillator_strength']['components']
    assert properties['oscillator_strength']['availability']
    for name in ['transition_dipole', 'transition_quadrupole']:
        for state in properties[name]['values']:
            assert not state['components'] and state['availability']


def assert_subworld_timing_scopes(attempts, groups, total_ranks):
    """Each participating World measures its own work; run totals are independent."""
    def walk(node):
        yield node
        for child in node.get('children', []):
            yield from walk(child)
    found = []
    for attempt in attempts:
        root = attempt['timing']
        nodes = [n for n in walk(root) if n.get('clock_scope') == 'subworld_driver_rank']
        for node in nodes:
            assert node['name'] == 'subworld_work'
            assert 0 <= node['subworld_id'] < groups
            assert node['mpi_ranks'] == total_ranks // groups
            assert node['threads_per_rank'] > 0
            assert 0 <= node['seconds'] <= root['seconds']
        found.extend(nodes)
    assert {n['subworld_id'] for n in found} == set(range(groups))
