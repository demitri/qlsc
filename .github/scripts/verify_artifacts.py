#!/usr/bin/env python3
'''
Check that a directory of release artifacts is complete and consistent.

Run by the "verify artifacts" job of .github/workflows/wheels.yml over every
wheel and sdist produced by a run, and by hand during a release over the
downloaded artifacts (see the release procedure in CLAUDE.md):

    python .github/scripts/verify_artifacts.py release-dist

It answers the question no single build job can: did every platform actually
produce every wheel it was supposed to, do they all carry the same version, and
(on a tag build, via the RELEASE_TAG environment variable) is that the version
the tag claims? A build job that silently produced fewer wheels than expected -
a matrix entry skipped, an interpreter dropped by a cibuildwheel upgrade - still
goes green on its own; this does not.

Filenames are parsed with packaging.utils rather than a local regex, so the tags
are interpreted by the same code that pip uses. That matters for cases a
hand-rolled pattern gets wrong: a wheel tagged "linux_x86_64" (auditwheel never
ran, so it is NOT installable as manylinux and must not be counted as one), a
"py3-none-any" wheel (the extension silently missing from the build), or an ABI
tag that disagrees with the interpreter tag.

Exits non-zero listing every discrepancy found. It never "repairs" what it finds
and never passes with a warning.

Requires the "packaging" library (a dependency of both pip and twine).
'''

import os
import sys
import pathlib

from packaging.utils import (
	parse_wheel_filename,
	parse_sdist_filename,
	canonicalize_name,
	InvalidWheelFilename,
	InvalidSdistFilename,
)
from packaging.version import Version, InvalidVersion

PROJECT = "qlsc"

# The interpreters and platforms the wheel matrix is expected to cover. Keep in
# step with "build" in source/pyproject.toml and the matrix in wheels.yml: this
# list is the statement of intent that those two are checked against.
EXPECTED_PYTHONS = ["cp310", "cp311", "cp312", "cp313", "cp314"]
EXPECTED_PLATFORMS = [
	"linux-manylinux-x86_64",
	"linux-musllinux-x86_64",
	"linux-manylinux-aarch64",
	"linux-musllinux-aarch64",
	"macos-arm64",
	"macos-x86_64",
	"windows-amd64",
]

# Architectures named inside a platform tag. Anything not listed here is
# reported rather than guessed at.
LINUX_ARCHES = ["x86_64", "aarch64", "i686", "armv7l", "ppc64le", "s390x", "riscv64"]
MACOS_ARCHES = ["arm64", "x86_64", "universal2", "universal", "intel"]


def platform_family(platform_tag: str) -> str:
	'''
	Reduce one wheel platform tag to a coarse "family" name.

	The families that a release is expected to contain are listed in
	EXPECTED_PLATFORMS; this function also names the ones that are NOT expected
	(an unrepaired Linux wheel, a pure-Python wheel, 32-bit Windows) so that
	they are reported as unexpected artifacts instead of being mistaken for a
	wheel that belongs in the release.

	:param platform_tag: a single platform tag, e.g. "manylinux_2_28_x86_64"
	:returns: a family name such as "linux-manylinux-x86_64"
	:raises ValueError: if the tag cannot be classified at all
	'''
	tag = platform_tag.lower()

	if tag == "any":
		# A wheel with no platform at all: the C extension did not make it into
		# the build. Never a valid qlsc artifact.
		return "pure-python-any"

	if tag.startswith("musllinux"):
		libc = "musllinux"
	elif tag.startswith("manylinux"):
		libc = "manylinux"
	elif tag.startswith("linux_"):
		# A raw linux_* tag means auditwheel never repaired the wheel. pip will
		# not install it from PyPI, and PyPI will not even accept it - so it
		# must not be counted as satisfying the manylinux slot.
		libc = "unrepaired"
	else:
		libc = None

	if libc is not None:
		for arch in LINUX_ARCHES:
			if tag.endswith(arch):
				return f"linux-{libc}-{arch}"
		raise ValueError(f"unrecognised Linux architecture in platform tag {platform_tag!r}")

	if tag.startswith("macosx"):
		for arch in MACOS_ARCHES:
			if tag.endswith(arch):
				return f"macos-{arch}"
		raise ValueError(f"unrecognised macOS architecture in platform tag {platform_tag!r}")

	if tag == "win_amd64":
		return "windows-amd64"
	if tag == "win_arm64":
		return "windows-arm64"
	if tag == "win32":
		return "windows-win32"

	raise ValueError(f"unrecognised platform tag {platform_tag!r}")


def check(dist: pathlib.Path, declared_version: str, release_tag: str = ""):
	'''
	Check the artifacts in a directory against what a release should contain.

	:param dist: directory holding the wheels and the sdist (files only)
	:param declared_version: the version in source/qlsc/version.py
	:param release_tag: the git tag being released, e.g. "v1.1.0"; empty when
		this is not a tag build, in which case no tag comparison is made
	:returns: a (problems, summary) tuple - problems is a list of strings, empty
		when the artifact set is exactly as expected
	'''
	problems = []
	found = {}		# (python tag, platform family) -> filename
	sdists = []
	versions = set()

	try:
		expected_version = Version(declared_version)
	except InvalidVersion:
		return ([f"the version in source/qlsc/version.py is not a valid version: {declared_version!r}"], {})

	entries = sorted(dist.iterdir(), key=lambda p: p.name)
	files = [p for p in entries if p.is_file()]
	for entry in entries:
		if not entry.is_file():
			problems.append(f"{entry.name}: not a file - the artifact directory must contain only wheels and the sdist")

	if not files:
		return ([f"no artifacts at all in {dist}"], {})

	for path in files:
		name = path.name
		if name.endswith(".whl"):
			try:
				project, version, _build, tags = parse_wheel_filename(name)
			except InvalidWheelFilename as exc:
				problems.append(f"{name}: not a valid wheel filename ({exc})")
				continue
			versions.add(version)
			if canonicalize_name(project) != PROJECT:
				problems.append(f"{name}: is a wheel for {project!r}, not {PROJECT!r}")
				continue

			interpreters = {tag.interpreter for tag in tags}
			abis = {tag.abi for tag in tags}
			families = set()
			for tag in tags:
				try:
					families.add(platform_family(tag.platform))
				except ValueError as exc:
					problems.append(f"{name}: {exc}")
			if len(interpreters) != 1:
				problems.append(f"{name}: carries more than one interpreter tag ({sorted(interpreters)})")
				continue
			if not families:
				continue
			if len(families) != 1:
				problems.append(f"{name}: platform tags disagree about the platform ({sorted(families)})")
				continue

			interpreter = interpreters.pop()
			family = families.pop()
			# Every wheel built here is a version-specific CPython extension, so
			# the ABI tag equals the interpreter tag (true for CPython 3.8+).
			# An abi3 or mismatched tag means the build changed shape and the
			# expectations below need revisiting deliberately.
			if abis != {interpreter}:
				problems.append(f"{name}: ABI tag {sorted(abis)} does not match interpreter tag {interpreter!r}")

			key = (interpreter, family)
			if key in found:
				problems.append(f"two wheels for {interpreter} {family}: {found[key]} and {name}")
			else:
				found[key] = name

		elif name.endswith(".tar.gz"):
			try:
				project, version = parse_sdist_filename(name)
			except InvalidSdistFilename as exc:
				problems.append(f"{name}: not a valid sdist filename ({exc})")
				continue
			if canonicalize_name(project) != PROJECT:
				problems.append(f"{name}: is an sdist for {project!r}, not {PROJECT!r}")
				continue
			versions.add(version)
			sdists.append(name)

		else:
			problems.append(f"{name}: not a qlsc wheel or sdist - unexpected artifact")

	# every expected interpreter/platform combination present, and nothing extra
	for python in EXPECTED_PYTHONS:
		for platform in EXPECTED_PLATFORMS:
			if (python, platform) not in found:
				problems.append(f"missing wheel: {python} {platform}")
	for python, platform in sorted(found):
		if python not in EXPECTED_PYTHONS or platform not in EXPECTED_PLATFORMS:
			problems.append(f"unexpected wheel: {python} {platform} ({found[(python, platform)]})")

	if len(sdists) != 1:
		problems.append(f"expected exactly one sdist, found {len(sdists)}: {sdists or 'none'}")

	# one version across the whole set, and it is the version in the source tree
	if len(versions) > 1:
		problems.append(f"artifacts carry more than one version: {sorted(str(v) for v in versions)}")
	elif versions and versions != {expected_version}:
		problems.append(f"artifacts are version {versions.pop()}, but source/qlsc/version.py says {expected_version}")

	# on a tag build, the tag has to agree with the version being shipped.
	#
	# The prefix is checked character by character rather than stripped: "v" is
	# the only thing allowed in front of the version, exactly once. str.lstrip
	# would not do - it removes EVERY leading "v", so a mistyped "vv1.1.0" tag
	# (which still matches the workflow's "v*" trigger, so it really does reach
	# here) would read as 1.1.0 and pass. Nor is it enough to remove one "v" and
	# hand the rest to Version: PEP 440 itself permits a leading "v", so
	# Version("v1.1.0") parses happily and compares equal to 1.1.0.
	tag = release_tag.strip()
	if tag:
		if not tag.startswith("v") or not tag[1:2].isdigit():
			problems.append(f"tag {tag} is not a version tag: it must be 'v' followed directly by the version, as in v{expected_version}")
		else:
			try:
				tagged_version = Version(tag[1:])
			except InvalidVersion:
				problems.append(f"tag {tag} is not a version tag: {tag[1:]!r} is not a valid version")
			else:
				if tagged_version != expected_version:
					problems.append(f"tag {tag} does not match the package version {expected_version}")

	summary = {
		"wheels": len(found),
		"sdists": len(sdists),
		"version": str(expected_version),
		"tag": tag,
		"found": found,
	}
	return (problems, summary)


def package_version(repo_root: pathlib.Path) -> str:
	'''
	Read __version__ out of source/qlsc/version.py without importing the
	package (which would need the compiled extension and NumPy).
	'''
	import re
	text = (repo_root / "source" / "qlsc" / "version.py").read_text()
	match = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", text)
	if match is None:
		raise SystemExit("could not find __version__ in source/qlsc/version.py")
	return match.group(1)


def main(argv):
	if len(argv) != 2:
		raise SystemExit(f"usage: {argv[0]} <directory of artifacts>")
	dist = pathlib.Path(argv[1])
	if not dist.is_dir():
		raise SystemExit(f"not a directory: {dist}")
	repo_root = pathlib.Path(__file__).resolve().parent.parent.parent

	release_tag = os.environ.get("RELEASE_TAG", "")
	problems, summary = check(dist, package_version(repo_root), release_tag)

	if summary:
		if summary["tag"]:
			print(f"tag {summary['tag']} checked against the package version {summary['version']}")
		else:
			print("not a tag build: skipping the tag/version comparison")
		print(f"{summary['wheels']} wheels and {summary['sdists']} sdist checked, version {summary['version']}")
		for python in EXPECTED_PYTHONS:
			built = [p for p in EXPECTED_PLATFORMS if (python, p) in summary["found"]]
			print(f"  {python}: {len(built)}/{len(EXPECTED_PLATFORMS)} platforms")

	if problems:
		print(f"\n{len(problems)} problem(s) with the artifact set:", file=sys.stderr)
		for problem in problems:
			print(f"  - {problem}", file=sys.stderr)
		return 1
	print("\nartifact set is complete and consistent")
	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))
