import sys
from pathlib import Path

import h5py
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from compare_madqc_hdf5 import classify_hdf5


def test_classifies_native_madness_archive(tmp_path):
    path = tmp_path / 'native.h5'
    with h5py.File(path, 'w') as handle:
        handle.create_dataset('madness_archive', data=[1, 2, 3])
        handle.create_group('meta')
    assert classify_hdf5(path) == 'native_madness_archive'


def test_classifies_structured_pymra_export(tmp_path):
    path = tmp_path / 'orbital.mad.h5'
    with h5py.File(path, 'w') as handle:
        handle.create_dataset('coeffs', data=[1.0])
        handle.create_dataset('keys', data=[1])
        handle.create_group('meta')
    assert classify_hdf5(path) == 'pymra_structured'


def test_missing_file_is_unavailable(tmp_path):
    with pytest.raises(ValueError, match='does not exist'):
        classify_hdf5(tmp_path / 'missing.h5')
