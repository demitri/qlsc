
import os
import sqlite3
import stat

import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex

# The "key" column carries an index of its own. Without it, every lookup or removal
# by key scans the whole table, which is the difference between a constant-time
# operation and one proportional to the size of the catalog.

TABLE = "qlsc_ipix"


def _make_index_file(path, n=500):
	index = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
	index.add_points(ra=np.linspace(0., 300., n),
					 dec=np.linspace(-80., 80., n),
					 keys=[f"k{i}" for i in range(n)])
	index.close()
	return path


def _index_names(path):
	connection = sqlite3.connect(path)
	try:
		return {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='index'")}
	finally:
		connection.close()


def _query_plan(path, sql, parameters):
	connection = sqlite3.connect(path)
	try:
		return " ".join(str(row) for row in connection.execute("EXPLAIN QUERY PLAN " + sql, parameters))
	finally:
		connection.close()


def _drop_key_index(path):
	'''
	Turn a current index file into one written before the key index existed.
	'''
	connection = sqlite3.connect(path)
	try:
		connection.execute("DROP INDEX IF EXISTS key_idx")
		connection.commit()
	finally:
		connection.close()


def test_new_index_file_has_a_key_index(tmp_path):
	'''
	Test that a newly created index file carries an index on the key column.
	'''
	path = _make_index_file(tmp_path / "new.qlscidx")
	assert "key_idx" in _index_names(path)


def test_key_queries_use_the_index(tmp_path):
	'''
	Test the property the index exists for: the queries behind "in", find_by_key
	and remove_point(key=...) search the index rather than scanning the table.

	This is the mechanical guard - a future change to those queries that stops them
	being able to use the index shows up here rather than as an unexplained
	slowdown on a large catalog.
	'''
	path = _make_index_file(tmp_path / "plans.qlscidx")

	statements = [
		(f"SELECT 1 FROM {TABLE} WHERE key=? LIMIT 1", ("k1",)),
		(f"SELECT ra,dec FROM {TABLE} WHERE key=?", ("k1",)),
		(f"DELETE FROM {TABLE} WHERE key=?", ("k1",))
	]

	for sql, parameters in statements:
		plan = _query_plan(path, sql, parameters)
		assert "key_idx" in plan, f"query does not use the key index: {sql} -> {plan}"
		assert "SCAN" not in plan, f"query still scans the table: {sql} -> {plan}"


def test_old_index_file_is_upgraded_on_open(tmp_path):
	'''
	Test the migration path: an index file written before the key index existed
	(qlsc 1.1 and earlier) gains it when it is opened, and its data is untouched.
	'''
	path = _make_index_file(tmp_path / "old.qlscidx")
	_drop_key_index(path)
	assert "key_idx" not in _index_names(path), "the fixture did not produce a pre-change file"

	index = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
	try:
		assert "key_idx" in _index_names(path)
		assert len(index) == 500
		assert "k42" in index
		assert index.find_by_key("k42")["ra"][0] == pytest.approx(np.linspace(0., 300., 500)[42])
	finally:
		index.close()

	# and the upgrade is persistent, not an artifact of the open connection
	assert "key_idx" in _index_names(path)


def test_old_index_file_works_before_and_after_the_upgrade(tmp_path):
	'''
	Test that the answers do not depend on the index being there: the same lookups
	and removals against a file with the key index dropped return the same results.
	'''
	path = _make_index_file(tmp_path / "answers.qlscidx", n=50)
	_drop_key_index(path)

	index = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
	try:
		assert sorted(index.keys()) == sorted(f"k{i}" for i in range(50))
		assert index.remove_point(key="k7") == 1
		assert "k7" not in index
		assert len(index) == 49
	finally:
		index.close()


def test_read_only_index_file_still_opens(tmp_path):
	'''
	Test that an index file that cannot be written to still opens, with a warning,
	rather than failing because the index could not be added.

	The index is an optimization; refusing to open a file that opened fine in the
	previous release would be the worse failure. The warning is what keeps it from
	being a silent difference.
	'''
	path = _make_index_file(tmp_path / "readonly.qlscidx", n=20)
	_drop_key_index(path)
	os.chmod(path, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
	if os.access(path, os.W_OK):
		# root (cibuildwheel's manylinux containers run the tests as root) can write a
		# file whatever its mode, so the read-only file this test is about cannot be
		# made here; the run on any other account covers it
		os.chmod(path, stat.S_IRUSR | stat.S_IWUSR)
		pytest.skip("cannot make a read-only file for this user (running as root?)")

	try:
		with _capture_warning() as captured:
			index = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
		try:
			assert len(index) == 20
			assert "k3" in index
			assert index.find_by_key("k3") is not None
		finally:
			index.close()

		assert "key_idx" not in _index_names(path), "nothing should have been written"
		assert captured, "opening a read-only index file must say why the key index is missing"
		assert "key" in captured[0].lower()
	finally:
		os.chmod(path, stat.S_IRUSR | stat.S_IWUSR)


class _capture_warning:
	'''
	Collect the messages qlsc logs at warning level while inside the block.
	'''
	def __init__(self):
		import logging
		self.records = list()
		self.logger = logging.getLogger("qlsc")

		captor = self

		class _Handler(logging.Handler):
			def emit(self, record):
				captor.records.append(record.getMessage())

		self.handler = _Handler(level=logging.WARNING)

	def __enter__(self):
		self.logger.addHandler(self.handler)
		self.previous_level = self.logger.level
		self.logger.setLevel(10) # DEBUG, so the warning is not filtered out
		return self.records

	def __exit__(self, exc_type, exc_value, traceback):
		self.logger.removeHandler(self.handler)
		self.logger.setLevel(self.previous_level)
		return False
