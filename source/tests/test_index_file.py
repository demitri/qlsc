
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex

def test_file_index_reopen(tmp_path):
	'''
	Test creating a file-based index, closing it, and reopening it.
	'''
	path = tmp_path / "test.qlscidx"

	idx = QLSCIndex(qlsc=QLSC(depth=30), filepath=path, name="test index")
	idx.add_point(12., 34.)
	idx.add_point(56., -78.)
	idx.close()

	idx2 = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
	assert idx2.number_of_points == 2
	assert idx2.name == "test index"
	matches = idx2.radial_query(ra=12., dec=34., radius=0.1)
	assert matches.shape == (2,)

def test_file_index_reopen_depth_mismatch(tmp_path):
	'''
	Test reopening an existing index file with a QLSC at a different depth:
	the index's stored depth wins.

	Regression test: this path raised NameError (undefined variable) when
	replacing the provided QLSC with one at the stored depth.
	'''
	path = tmp_path / "test.qlscidx"

	idx = QLSCIndex(qlsc=QLSC(depth=25), filepath=path)
	idx.add_point(12., 34.)
	idx.close()

	idx2 = QLSCIndex(qlsc=QLSC(depth=30), filepath=path) # depth differs from the file's
	assert idx2.qlsc.depth == 25
	assert idx2.number_of_points == 1
	matches = idx2.radial_query(ra=12., dec=34., radius=0.1)
	assert matches.shape == (2,)
