'''
Tests for verify_artifacts.py, the release gate.

The gate is only worth having if it actually fails on a broken artifact set, so
each test builds a synthetic set (empty files - only the names are read) and
asserts on the specific problem reported. A gate that silently passes everything
would look identical in CI to a gate that works; these tests are what tells the
two apart.

Run by the "verify script tests" job of .github/workflows/wheels.yml, and by
hand with: pytest .github/scripts
'''

import pathlib
import pytest

import verify_artifacts
from verify_artifacts import check, platform_family, EXPECTED_PYTHONS, EXPECTED_PLATFORMS

VERSION = "1.1.0"

# one representative filename per expected platform family
PLATFORM_TAGS = {
	"linux-manylinux-x86_64": "manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64",
	"linux-musllinux-x86_64": "musllinux_1_2_x86_64",
	"linux-manylinux-aarch64": "manylinux2014_aarch64.manylinux_2_17_aarch64.manylinux_2_28_aarch64",
	"linux-musllinux-aarch64": "musllinux_1_2_aarch64",
	"macos-arm64": "macosx_11_0_arm64",
	"macos-x86_64": "macosx_10_9_x86_64",
	"windows-amd64": "win_amd64",
}


def wheel_name(python: str, family: str, version: str = VERSION) -> str:
	return f"qlsc-{version}-{python}-{python}-{PLATFORM_TAGS[family]}.whl"


def make_set(tmp_path: pathlib.Path, version: str = VERSION, omit=(), extra=()) -> pathlib.Path:
	'''
	Build a complete, valid artifact set, minus anything in "omit" and plus any
	literal filenames in "extra".
	'''
	dist = tmp_path / "dist"
	dist.mkdir(exist_ok=True)
	for python in EXPECTED_PYTHONS:
		for family in EXPECTED_PLATFORMS:
			if (python, family) in omit:
				continue
			(dist / wheel_name(python, family, version)).write_bytes(b"")
	if "sdist" not in omit:
		(dist / f"qlsc-{version}.tar.gz").write_bytes(b"")
	for name in extra:
		(dist / name).write_bytes(b"")
	return dist


def test_complete_set_has_no_problems(tmp_path):
	problems, summary = check(make_set(tmp_path), VERSION)
	assert problems == []
	assert summary["wheels"] == len(EXPECTED_PYTHONS) * len(EXPECTED_PLATFORMS) == 35
	assert summary["sdists"] == 1


def test_missing_wheel_is_reported(tmp_path):
	dist = make_set(tmp_path, omit=[("cp312", "windows-amd64")])
	problems, _ = check(dist, VERSION)
	assert problems == ["missing wheel: cp312 windows-amd64"]


def test_missing_sdist_is_reported(tmp_path):
	problems, _ = check(make_set(tmp_path, omit=["sdist"]), VERSION)
	assert any("expected exactly one sdist" in p for p in problems)


def test_two_sdists_are_reported(tmp_path):
	# both names are valid sdist filenames, so this exercises the duplicate
	# branch and not filename validation
	dist = make_set(tmp_path, extra=["qlsc-1.1.1.tar.gz"])
	problems, _ = check(dist, VERSION)
	assert any(p.startswith("expected exactly one sdist, found 2") for p in problems)


def test_malformed_sdist_filename_is_reported(tmp_path):
	dist = make_set(tmp_path, extra=["qlsc-1.1.0.tar.gz.tar.gz"])
	problems, _ = check(dist, VERSION)
	assert any("not a valid sdist filename" in p for p in problems)


def test_duplicate_wheel_for_one_slot_is_reported(tmp_path):
	# a second cp310 manylinux x86_64 wheel, tagged with a different (but still
	# manylinux x86_64) platform string
	dist = make_set(tmp_path, extra=["qlsc-1.1.0-cp310-cp310-manylinux_2_17_x86_64.whl"])
	problems, _ = check(dist, VERSION)
	assert any(p.startswith("two wheels for cp310 linux-manylinux-x86_64") for p in problems)


def test_malformed_wheel_filename_is_reported(tmp_path):
	dist = make_set(tmp_path, extra=["qlsc-1.1.0-cp313.whl"])
	problems, _ = check(dist, VERSION)
	assert any("not a valid wheel filename" in p for p in problems)


def test_stray_file_is_reported(tmp_path):
	dist = make_set(tmp_path, extra=["build-notes.txt"])
	problems, _ = check(dist, VERSION)
	assert any("unexpected artifact" in p for p in problems)


def test_subdirectory_is_reported(tmp_path):
	# the failure mode of a release procedure that flattens "gh run download"
	# output badly: the per-artifact directories are left behind
	dist = make_set(tmp_path)
	(dist / "wheels-linux-x86_64").mkdir()
	problems, _ = check(dist, VERSION)
	assert any("not a file" in p for p in problems)


def test_unrepaired_linux_wheel_does_not_satisfy_the_manylinux_slot(tmp_path):
	'''
	A linux_x86_64 wheel is what comes out of a build auditwheel never repaired.
	PyPI rejects it and pip will not install it, so it must be reported rather
	than counted as the manylinux wheel for that interpreter.
	'''
	dist = make_set(tmp_path, omit=[("cp311", "linux-manylinux-x86_64")],
	                extra=["qlsc-1.1.0-cp311-cp311-linux_x86_64.whl"])
	problems, _ = check(dist, VERSION)
	assert "missing wheel: cp311 linux-manylinux-x86_64" in problems
	assert any("unexpected wheel: cp311 linux-unrepaired-x86_64" in p for p in problems)


def test_pure_python_wheel_is_reported(tmp_path):
	'''A py3-none-any wheel means the C extension never got built in.'''
	dist = make_set(tmp_path, extra=["qlsc-1.1.0-py3-none-any.whl"])
	problems, _ = check(dist, VERSION)
	assert any("pure-python-any" in p for p in problems)


def test_abi_tag_mismatch_is_reported(tmp_path):
	dist = make_set(tmp_path, extra=["qlsc-1.1.0-cp313-cp312-win_arm64.whl"])
	problems, _ = check(dist, VERSION)
	assert any("ABI tag" in p for p in problems)


def test_wheel_for_another_project_is_reported(tmp_path):
	dist = make_set(tmp_path, extra=["other-1.0.0-cp310-cp310-win_amd64.whl"])
	problems, _ = check(dist, VERSION)
	assert any("not 'qlsc'" in p for p in problems)


def test_version_disagreeing_with_version_py_is_reported(tmp_path):
	problems, _ = check(make_set(tmp_path, version="1.2.0"), VERSION)
	assert any("but source/qlsc/version.py says 1.1.0" in p for p in problems)


def test_mixed_versions_are_reported(tmp_path):
	dist = make_set(tmp_path, extra=["qlsc-1.2.0-cp310-cp310-win_arm64.whl"])
	problems, _ = check(dist, VERSION)
	assert any("more than one version" in p for p in problems)


def test_matching_tag_passes(tmp_path):
	problems, summary = check(make_set(tmp_path), VERSION, release_tag="v1.1.0")
	assert problems == []
	assert summary["tag"] == "v1.1.0"


def test_mismatched_tag_is_reported(tmp_path):
	problems, _ = check(make_set(tmp_path), VERSION, release_tag="v9.9.9")
	assert problems == ["tag v9.9.9 does not match the package version 1.1.0"]


def test_non_version_tag_is_reported(tmp_path):
	problems, _ = check(make_set(tmp_path), VERSION, release_tag="release-candidate")
	assert any("not a version tag" in p for p in problems)


@pytest.mark.parametrize("tag", ["vv1.1.0", "vvv1.1.0", "1.1.0", "version1.1.0", "v-1.1.0", "v"])
def test_malformed_version_tags_are_reported(tag, tmp_path):
	'''
	"v" prefixed exactly once, or it is not a release tag.

	"vv1.1.0" is the case worth spelling out: it still matches the workflow's
	"v*" push trigger, so it reaches this check for real, and two plausible
	implementations both wave it through - str.lstrip("v") strips both letters,
	and PEP 440 lets Version() accept a leading "v" of its own, so peeling one
	off and parsing the rest also yields 1.1.0.
	'''
	problems, _ = check(make_set(tmp_path), VERSION, release_tag=tag)
	assert any("not a version tag" in p for p in problems), f"{tag!r} was accepted as a release tag"


def test_empty_directory_is_reported(tmp_path):
	dist = tmp_path / "empty"
	dist.mkdir()
	problems, _ = check(dist, VERSION)
	assert any("no artifacts at all" in p for p in problems)


@pytest.mark.parametrize("tag,expected", [
	("manylinux_2_28_x86_64", "linux-manylinux-x86_64"),
	("manylinux2014_aarch64", "linux-manylinux-aarch64"),
	("musllinux_1_2_x86_64", "linux-musllinux-x86_64"),
	("linux_x86_64", "linux-unrepaired-x86_64"),
	("macosx_11_0_arm64", "macos-arm64"),
	("macosx_10_13_x86_64", "macos-x86_64"),
	("macosx_10_9_universal2", "macos-universal2"),
	("win_amd64", "windows-amd64"),
	("win32", "windows-win32"),
	("any", "pure-python-any"),
])
def test_platform_family(tag, expected):
	assert platform_family(tag) == expected


def test_platform_family_rejects_the_unclassifiable():
	with pytest.raises(ValueError):
		platform_family("manylinux_2_28_sparc")
	with pytest.raises(ValueError):
		platform_family("plan9_386")


def repo_root() -> pathlib.Path:
	return pathlib.Path(verify_artifacts.__file__).resolve().parent.parent.parent


def test_expected_pythons_match_the_cibuildwheel_configuration():
	'''
	The expected interpreter list is a duplicate of "build" in
	source/pyproject.toml; if the two drift apart the gate checks for wheels
	nobody builds, or waves through a dropped interpreter.
	'''
	pyproject = (repo_root() / "source" / "pyproject.toml").read_text()
	for python in EXPECTED_PYTHONS:
		assert f'"{python}-*"' in pyproject, f"{python} is expected here but not built by cibuildwheel"
	assert pyproject.count("cp3") >= len(EXPECTED_PYTHONS)


def test_expected_platforms_match_the_workflow_matrix():
	'''
	The same drift check for platforms: EXPECTED_PLATFORMS duplicates what the
	wheels.yml build matrix actually covers. Deleting a matrix entry without
	touching this list would leave the gate demanding wheels nobody builds
	(every run then red for the wrong reason); adding one without touching the
	list would let a whole platform go unchecked.

	Requires PyYAML, so that the matrix is read the way Actions reads it.
	'''
	import yaml
	workflow = yaml.safe_load((repo_root() / ".github" / "workflows" / "wheels.yml").read_text())
	matrix = workflow["jobs"]["wheels"]["strategy"]["matrix"]["include"]

	covered = set()
	for entry in matrix:
		os_label, arch = entry["os"], entry["archs"]
		if os_label.startswith("ubuntu"):
			# each Linux runner produces both a manylinux and a musllinux wheel
			covered.add(f"linux-manylinux-{arch}")
			covered.add(f"linux-musllinux-{arch}")
		elif os_label.startswith("macos"):
			covered.add(f"macos-{arch}")
		elif os_label.startswith("windows"):
			covered.add(f"windows-{arch.lower()}")
		else:
			raise AssertionError(f"unrecognised runner {os_label!r} in the wheels matrix")

	assert covered == set(EXPECTED_PLATFORMS), (
		f"the wheels.yml matrix and EXPECTED_PLATFORMS disagree: "
		f"built but not expected {sorted(covered - set(EXPECTED_PLATFORMS))}, "
		f"expected but not built {sorted(set(EXPECTED_PLATFORMS) - covered)}")
