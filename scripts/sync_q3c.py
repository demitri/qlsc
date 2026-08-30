#!/usr/bin/env python3

'''
Sync the vendored Q3C C sources in this repository to an upstream Q3C release.

This mechanizes the procedure documented in the "Vendored Q3C C code" section of
CLAUDE.md. It fetches a tag from https://github.com/segasai/q3c, copies the four
vendored files verbatim, re-applies whichever of the local patches upstream still
lacks, updates the Q3C_VERSION string in setup.py, checks that the C interface the
wrapper depends on has not changed, and runs the test suite.

Usage
-----

	# see what syncing to q3c v2.0.6 would do, changing nothing
	python3 scripts/sync_q3c.py v2.0.6 --dry-run

	# do it
	python3 scripts/sync_q3c.py v2.0.6

	# re-run the current version against itself as an end-to-end check
	# (expected result: "no changes", both patches still needed)
	python3 scripts/sync_q3c.py v2.0.5 --dry-run

	# work from already-unpacked upstream trees instead of downloading
	python3 scripts/sync_q3c.py v2.0.6 --upstream-dir /tmp/q3c-2.0.6 \
	                                   --previous-dir /tmp/q3c-2.0.5

Options: --dry-run (report only), --from-tag (the upstream tag the vendored files
currently come from; defaults to the Q3C_VERSION in setup.py), --upstream-dir /
--previous-dir (use a local tree instead of downloading), --skip-tests (do not
reinstall and run pytest).

What it does, in order
----------------------

1. Fetches the target tag and the tag the tree currently sits on. Both are needed:
   the previous release is the baseline that proves the vendored files carry the
   two documented local patches and *nothing else*, and its common.h is the "before"
   side of the interface diff.
2. Confirms the vendored files are exactly previous-upstream plus the known patches.
   An undocumented local divergence stops the run: it would be silently overwritten.
3. For each patch, decides whether the new upstream release already contains the fix
   (drop the patch) or still carries the bug (re-apply it). The decision is made by
   exact string match on the surrounding code, and an ambiguous result - the anchor
   text found the wrong number of times, or neither the original nor the patched form
   present - is a hard error, never a fuzzy match. Fuzzing a patch into a stack buffer
   overflow fix is precisely the failure this script exists to prevent.
4. Diffs the C interface: every identifier the wrapper (qlsc_c_module.c) uses that is
   declared in the old common.h is looked up in the new one. A changed prototype,
   struct, or macro - or one that disappeared - stops the run.
5. Derives the Q3C_VERSION string from the tag, cross-checks it against upstream's
   own q3c.control, and updates setup.py.
6. Writes the four files, reinstalls the package (the test suite imports the INSTALLED
   package, not the working tree), and runs the guard tests followed by the full suite.

Anything this script cannot decide it refuses to guess at. After a successful run there
is still human work to do; it prints the checklist.
'''

import argparse
import difflib
import enum
import pathlib
import re
import subprocess
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCE_DIR = REPO_ROOT / "source"
VENDOR_DIR = SOURCE_DIR / "qlsc" / "c_code" / "q3c"
WRAPPER_PATH = SOURCE_DIR / "qlsc" / "c_code" / "qlsc_c_module.c"
SETUP_PATH = SOURCE_DIR / "setup.py"
TESTS_DIR = SOURCE_DIR / "tests"

UPSTREAM_REPO_URL = "https://github.com/segasai/q3c"
UPSTREAM_ISSUES_URL = "https://github.com/segasai/q3c/issues"

# The files vendored into source/qlsc/c_code/q3c/ . q3c.c (the PostgreSQL layer)
# is deliberately not among them.
VENDORED_FILES = ("common.h", "my_bits.h", "q3c_poly.c", "q3cube.c")

# The header carrying the interface the wrapper compiles against.
HEADER_FILE = "common.h"

# Macros used by the wrapper that are not defined in the vendored sources, so the
# interface diff must not expect to find them there.
EXTERNALLY_DEFINED_MACROS = frozenset({
	"Q3C_VERSION",	 # -DQ3C_VERSION, from setup.py's define_macros
	"Q3C_STANDALONE" # -DQ3C_STANDALONE, from setup.py's define_macros
})


class SyncError(Exception):
	'''
	A condition the sync cannot resolve on its own; the message says what a human must do.
	'''
	pass


class PatchState(enum.Enum):
	'''
	Whether a local patch is still needed against a given upstream source file.
	'''
	NEEDED = "needed"                      # upstream still carries the bug; apply the patch
	ALREADY_UPSTREAM = "already upstream"  # upstream now carries the fix; drop the patch


class Patch:
	'''
	One local divergence from upstream, expressed as an exact-match text substitution.

	Written out rather than declared with :mod:`dataclasses`, which arrived in Python
	3.7 while the package declares ``python_requires='>=3.6'``; the tests import this
	script, so a dataclass here would fail the suite on a supported interpreter.

	:param name: short identifier used in this script's output
	:param filename: the vendored file the patch applies to
	:param upstream_text: the text as upstream ships it (the "before" side)
	:param patched_text: the text as this repository ships it (the "after" side)
	:param occurrences: how many times the text appears in the file; any other count is an error
	:param issue: the upstream issue this patch works around
	:param guard_test: the test in source/tests/test_source_invariants.py that guards it
	'''
	__slots__ = ("name", "filename", "upstream_text", "patched_text", "occurrences", "issue", "guard_test")

	def __init__(self, name, filename, upstream_text, patched_text, occurrences, issue, guard_test):
		self.name = name
		self.filename = filename
		self.upstream_text = upstream_text
		self.patched_text = patched_text
		self.occurrences = occurrences
		self.issue = issue
		self.guard_test = guard_test

	def __repr__(self):
		return "Patch(name={!r}, filename={!r}, occurrences={!r})".format(self.name, self.filename, self.occurrences)


# ------------------------------------------------------------------
# The local patches. Both are documented in CLAUDE.md; keep the two in step.
# ------------------------------------------------------------------

POLY_POINTS_PATCH = Patch(
	name="q3c_poly.c points[8]",
	filename="q3c_poly.c",
	upstream_text="    q3c_coord_t points[4];\n",
	patched_text="    q3c_coord_t points[8];\n",
	occurrences=1,
	issue=f"{UPSTREAM_ISSUES_URL}/50",
	guard_test="test_q3c_poly_points_buffer_patch")

RES_DEPTH_PATCH = Patch(
	name="q3cube.c res_depth log2",
	filename="q3cube.c",
	upstream_text="        res_depth = nside / n0;\n",
	patched_text=(
		"        /* Local qlsc patch (upstream q3c issue #51): res_depth counts\n"
		"         * classification passes, of which only the first res_depth - 1 subdivide,\n"
		"         * so the safe value is log2(nside / n0) + 1, not the raw ratio. For tiny\n"
		"         * query regions the ratio exceeds the true remaining depth and the\n"
		"         * subdivision below overflows (2 * 2^30 in a signed int), producing\n"
		"         * garbage coordinates; and at least one classification pass must always\n"
		"         * run, or the squares' status fields are read uninitialized. Behavior is\n"
		"         * unchanged whenever the ratio >= 8, since both values are then capped\n"
		"         * to Q3C_MAX_DEPTH.\n"
		"         */\n"
		"        res_depth = 1;\n"
		"        for (q3c_ipix_t subdiv_ratio = nside / n0; subdiv_ratio > 1; subdiv_ratio >>= 1)\n"
		"        {\n"
		"            res_depth++;\n"
		"        }\n"),
	occurrences=3,
	issue=f"{UPSTREAM_ISSUES_URL}/51",
	guard_test="test_q3cube_res_depth_patch")

PATCHES = (POLY_POINTS_PATCH, RES_DEPTH_PATCH)


# ------------------------------------------------------------------
# Patch detection and application (pure text; tested in
# source/tests/test_sync_q3c.py)
# ------------------------------------------------------------------

def patch_state(patch:Patch, source:str) -> PatchState:
	'''
	Decide whether `patch` still needs to be applied to the given source text.

	The two accepted states are exact: the upstream text present the expected number
	of times and the patched text absent (the patch is still needed), or the patched
	text present the expected number of times and the upstream text absent (upstream
	has adopted the fix). Anything else - a changed count, a partial match, neither
	form present - raises, because there is no safe way to guess where the patch goes.

	:param patch: the patch to test for
	:param source: the contents of the file the patch applies to
	:returns: a :class:`PatchState`
	'''
	n_upstream = source.count(patch.upstream_text)
	n_patched = source.count(patch.patched_text)

	if n_patched == patch.occurrences and n_upstream == 0:
		return PatchState.ALREADY_UPSTREAM
	if n_upstream == patch.occurrences and n_patched == 0:
		return PatchState.NEEDED

	raise SyncError(
		f"The local patch '{patch.name}' cannot be placed in '{patch.filename}'.\n"
		f"Expected either {patch.occurrences} occurrence(s) of the upstream text and none of the\n"
		f"patched text, or the reverse; found {n_upstream} upstream and {n_patched} patched.\n"
		f"The surrounding code has changed, so this script will not guess where the patch belongs.\n"
		f"WHAT TO DO: read {patch.issue} and the new '{patch.filename}', decide by hand whether the\n"
		f"bug is still present, re-apply or drop the patch, and update both CLAUDE.md and the guard\n"
		f"test '{patch.guard_test}' in source/tests/test_source_invariants.py.")


def apply_patch(patch:Patch, source:str) -> tuple:
	'''
	Apply `patch` to the given source text if it is still needed.

	:param patch: the patch to apply
	:param source: the contents of the file the patch applies to
	:returns: a tuple of (new source text, :class:`PatchState`)
	'''
	state = patch_state(patch, source)
	if state is PatchState.ALREADY_UPSTREAM:
		return source, state
	return source.replace(patch.upstream_text, patch.patched_text), state


def unapply_patch(patch:Patch, source:str) -> str:
	'''
	Reverse `patch`, returning the text as upstream ships it.

	Used to check the vendored files against the upstream release they came from:
	reversing the known patches must leave the upstream file byte for byte.

	:param patch: the patch to reverse
	:param source: the contents of the patched (vendored) file
	:returns: the source text with the patch removed
	'''
	state = patch_state(patch, source)
	if state is PatchState.NEEDED:
		# the patch is not in this text at all
		raise SyncError(
			f"The local patch '{patch.name}' is not present in the vendored '{patch.filename}',\n"
			f"so it cannot be reversed. The vendored file is not the patched file this script expects.")
	return source.replace(patch.patched_text, patch.upstream_text)


def patches_for(filename:str) -> tuple:
	'''
	Return the patches that apply to the named vendored file.

	:param filename: the base name of a vendored file
	:returns: a tuple of :class:`Patch`
	'''
	return tuple(p for p in PATCHES if p.filename == filename)


# ------------------------------------------------------------------
# C source parsing for the interface diff
# ------------------------------------------------------------------

def strip_c_comments(text:str) -> str:
	'''
	Remove C comments from source text, leaving string and character literals intact.

	A regex would corrupt any string literal containing "//" (the wrapper has one: a
	URL in an error message), silently dropping identifiers that follow it on the same
	line. Identifiers dropped without notice are exactly what the interface diff must
	never do, so this walks the text instead.

	:param text: C source text
	:returns: the same text with comments replaced by a single space
	'''
	out = []
	i = 0
	n = len(text)
	while i < n:
		c = text[i]
		if c == '"' or c == "'":
			quote = c
			out.append(c)
			i += 1
			while i < n:
				if text[i] == "\\":
					out.append(text[i:i+2])
					i += 2
					continue
				out.append(text[i])
				i += 1
				if text[i-1] == quote:
					break
			continue
		if c == "/" and i + 1 < n and text[i+1] == "/":
			while i < n and text[i] != "\n":
				i += 1
			out.append(" ")
			continue
		if c == "/" and i + 1 < n and text[i+1] == "*":
			end = text.find("*/", i + 2)
			if end == -1:
				raise SyncError("Unterminated C block comment; the source file is not valid C.")
			# a space so the tokens either side of the comment do not run together,
			# and the newlines so the line-oriented handling below still lines up
			out.append(" " + "\n" * text.count("\n", i, end))
			i = end + 2
			continue
		out.append(c)
		i += 1
	return "".join(out)


def _join_continued_lines(lines:list) -> list:
	'''
	Join lines ending in a backslash with the line that follows, as the preprocessor does.

	:param lines: source lines without their terminating newlines
	:returns: the joined lines
	'''
	joined = []
	buffer = None
	for line in lines:
		part = line[:-1] if line.endswith("\\") else line
		buffer = part if buffer is None else buffer + " " + part
		if not line.endswith("\\"):
			joined.append(buffer)
			buffer = None
	if buffer is not None:
		joined.append(buffer)
	return joined


def _normalize(text:str) -> str:
	'''
	Collapse all runs of whitespace to a single space so that reformatting alone
	(upstream reruns clang-format now and then) is not reported as an interface change.

	:param text: any text
	:returns: the whitespace-normalized text
	'''
	return " ".join(text.split())


def _split_statements(text:str) -> list:
	'''
	Split C source into statements at semicolons that sit outside any braces, so that
	a struct or enum body stays with its declaration.

	:param text: comment- and directive-free C source text
	:returns: a list of statement strings
	'''
	statements = []
	depth = 0
	start = 0
	for i, c in enumerate(text):
		if c == "{":
			depth += 1
		elif c == "}":
			depth -= 1
		elif c == ";" and depth == 0:
			statements.append(text[start:i])
			start = i + 1
	statements.append(text[start:])
	return statements


_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")
_DEFINE_RE = re.compile(r"^\s*#\s*define\s+([A-Za-z_]\w*)")
_DIRECTIVE_RE = re.compile(r"^\s*#")
_STRUCT_DEF_RE = re.compile(r"\bstruct\s+([A-Za-z_]\w*)\s*\{")


def interface_index(header_text:str) -> dict:
	'''
	Index everything a header defines, by name.

	Covers macros, typedefs, struct/enum definitions and function declarations. A name
	can appear more than once (the q3c header defines its math macros twice, once in
	each branch of an #ifdef), so every name maps to the list of its definitions in
	source order; comparing the lists catches a change in either branch.

	:param header_text: the contents of a C header
	:returns: a dict mapping each defined name to a list of normalized definition strings
	'''
	index = {}

	def record(name, definition):
		index.setdefault(name, []).append(_normalize(definition))

	text = strip_c_comments(header_text)

	# macros first, then drop every preprocessor line so the remaining text is plain C
	body_lines = []
	for line in _join_continued_lines(text.split("\n")):
		match = _DEFINE_RE.match(line)
		if match:
			record(match.group(1), line)
			continue
		if _DIRECTIVE_RE.match(line):
			continue
		body_lines.append(line)

	for statement in _split_statements("\n".join(body_lines)):
		normalized = _normalize(statement)
		if not normalized:
			continue
		struct_def = _STRUCT_DEF_RE.search(normalized)
		if struct_def and normalized.rstrip().endswith("}"):
			# "struct name { ... }" - a definition, not a "struct name *" parameter
			record(struct_def.group(1), normalized)
			continue
		# strip any brace block so the declared name is found in the declarator, not the body
		declarator = re.sub(r"\{.*\}", " ", normalized)
		paren = declarator.find("(")
		if paren != -1:
			names = _IDENTIFIER_RE.findall(declarator[:paren])
		else:
			names = _IDENTIFIER_RE.findall(declarator)
		if not names:
			continue
		record(names[-1], normalized)

	return index


def wrapper_identifiers(wrapper_text:str) -> set:
	'''
	Every identifier the wrapper source mentions outside comments and string literals.

	:param wrapper_text: the contents of qlsc_c_module.c
	:returns: the set of identifiers
	'''
	stripped = strip_c_comments(wrapper_text)
	# remove string and character literals: an identifier inside a message is not a use
	stripped = re.sub(r'"(?:[^"\\]|\\.)*"', " ", stripped)
	stripped = re.sub(r"'(?:[^'\\]|\\.)*'", " ", stripped)
	return set(_IDENTIFIER_RE.findall(stripped))


def wrapper_locally_defined(wrapper_text:str) -> set:
	'''
	The names the wrapper defines for itself with #define; these are not part of the
	interface it inherits from the vendored header.

	:param wrapper_text: the contents of qlsc_c_module.c
	:returns: the set of locally defined macro names
	'''
	names = set()
	for line in _join_continued_lines(strip_c_comments(wrapper_text).split("\n")):
		match = _DEFINE_RE.match(line)
		if match:
			names.add(match.group(1))
	return names


def interface_differences(old_header:str, new_header:str, wrapper_text:str) -> list:
	'''
	Compare the parts of the C interface the wrapper actually uses between two versions
	of the header.

	Only names the wrapper mentions AND the old header defines are compared, so upstream
	is free to change anything qlsc does not touch. A name that the old header defines
	and the new one does not is reported as removed.

	:param old_header: contents of common.h from the release the vendored files come from
	:param new_header: contents of common.h from the release being synced to
	:param wrapper_text: the contents of qlsc_c_module.c
	:returns: a list of human-readable difference descriptions; empty if nothing the wrapper uses changed
	'''
	old_index = interface_index(old_header)
	new_index = interface_index(new_header)
	used = wrapper_identifiers(wrapper_text)

	differences = []
	for name in sorted(used & set(old_index)):
		old_definitions = old_index[name]
		new_definitions = new_index.get(name)
		if new_definitions is None:
			differences.append(f"'{name}' is used by the wrapper but no longer defined in {HEADER_FILE}:\n"
							   + "\n".join(f"    was: {d}" for d in old_definitions))
		elif new_definitions != old_definitions:
			lines = [f"'{name}' changed in {HEADER_FILE}:"]
			lines += [f"    was: {d}" for d in old_definitions]
			lines += [f"    now: {d}" for d in new_definitions]
			differences.append("\n".join(lines))

	# Names the wrapper uses that the old header never defined are not part of this
	# interface at all: they come from Python, NumPy, the C library, or setup.py's
	# define_macros (EXTERNALLY_DEFINED_MACROS names the latter, for the reader).
	return differences


# ------------------------------------------------------------------
# Version handling
# ------------------------------------------------------------------

_TAG_RE = re.compile(r"^v?(\d+(?:\.\d+)*)$")
_SETUP_VERSION_RE = re.compile(r"""\(\s*'Q3C_VERSION'\s*,\s*'"(?P<version>[^"]*)"'\s*\)""")
_CONTROL_VERSION_RE = re.compile(r"^\s*default_version\s*=\s*'(?P<version>[^']*)'", re.MULTILINE)


def version_from_tag(tag:str) -> str:
	'''
	Derive the Q3C_VERSION string from an upstream tag name, e.g. "v2.0.5" -> "2.0.5".

	:param tag: an upstream git tag
	:returns: the version string
	'''
	match = _TAG_RE.match(tag.strip())
	if match is None:
		raise SyncError(f"The tag '{tag}' is not of the form 'v2.0.5'; the Q3C_VERSION string cannot be "
						f"derived from it mechanically. Set it in setup.py by hand and pass the exact tag name.")
	return match.group(1)


def setup_q3c_version(setup_text:str) -> str:
	'''
	Read the Q3C_VERSION string compiled into the extension from setup.py.

	:param setup_text: the contents of setup.py
	:returns: the version string
	'''
	match = _SETUP_VERSION_RE.search(setup_text)
	if match is None:
		raise SyncError("Could not find the ('Q3C_VERSION', '\"x.y.z\"') entry in setup.py's define_macros; "
						"it has been reworded, so this script cannot update it. Update it by hand.")
	return match.group("version")


def set_setup_q3c_version(setup_text:str, version:str) -> str:
	'''
	Return setup.py's text with the Q3C_VERSION string replaced.

	:param setup_text: the contents of setup.py
	:param version: the version string to write
	:returns: the updated text
	'''
	setup_q3c_version(setup_text) # raises if the entry is not there
	return _SETUP_VERSION_RE.sub(f"""('Q3C_VERSION', '"{version}"')""", setup_text, count=1)


def control_version(control_text:str) -> str:
	'''
	Read upstream's own version string from its q3c.control file.

	:param control_text: the contents of q3c.control
	:returns: the version string
	'''
	match = _CONTROL_VERSION_RE.search(control_text)
	if match is None:
		raise SyncError("Upstream's q3c.control has no default_version line; the tag's version string "
						"cannot be cross-checked.")
	return match.group("version")


# ------------------------------------------------------------------
# Fetching upstream
# ------------------------------------------------------------------

def download_upstream(tag:str, destination:pathlib.Path) -> pathlib.Path:
	'''
	Download and unpack an upstream q3c release tarball.

	:param tag: the upstream git tag, e.g. "v2.0.5"
	:param destination: a directory to unpack into
	:returns: the path of the unpacked source tree
	'''
	if "/" in tag or tag in ("", ".", ".."):
		raise SyncError(f"'{tag}' is not a usable git tag name.")
	url = f"{UPSTREAM_REPO_URL}/archive/refs/tags/{tag}.tar.gz"
	print(f"  fetching {url}")
	tarball = destination / "release.tar.gz"
	try:
		with urllib.request.urlopen(url, timeout=60) as response:
			tarball.write_bytes(response.read())
	except urllib.error.HTTPError as e:
		if e.code == 404:
			raise SyncError(f"No such upstream tag '{tag}' ({url} returned 404). Check the tag name at "
							f"{UPSTREAM_REPO_URL}/tags .") from e
		raise SyncError(f"Could not download {url}: {e}") from e
	except urllib.error.URLError as e:
		raise SyncError(f"Could not download {url}: {e}. If this machine has no network access, unpack the "
						f"release by hand and pass --upstream-dir.") from e

	unpacked = destination / "unpacked"
	unpacked.mkdir()
	root = unpacked.resolve()
	# Unpacked member by member rather than with extractall(): only the file contents
	# are wanted, so nothing here can act on a symlink, a device node or an archived
	# file mode, and there is no version-dependent extraction filter to reason about.
	with tarfile.open(tarball) as archive:
		for member in archive.getmembers():
			if not (member.isfile() or member.isdir()):
				raise SyncError(f"The upstream tarball contains '{member.name}', which is neither a file nor a "
								f"directory. Refusing to unpack it.")
			target = (unpacked / member.name).resolve()
			if root not in target.parents:
				raise SyncError(f"The upstream tarball contains a path escaping the unpack directory: "
								f"'{member.name}'. Refusing to unpack it.")
			if member.isdir():
				target.mkdir(parents=True, exist_ok=True)
				continue
			target.parent.mkdir(parents=True, exist_ok=True)
			source = archive.extractfile(member)
			if source is None:
				raise SyncError(f"The upstream tarball's entry '{member.name}' has no readable contents.")
			with source, open(target, "wb") as out:
				out.write(source.read())

	roots = [p for p in unpacked.iterdir() if p.is_dir()]
	if len(roots) != 1:
		raise SyncError(f"Expected exactly one top level directory in the upstream tarball for '{tag}'; "
						f"found {len(roots)}.")
	return roots[0]


def upstream_tree(tag:str, override:str, workdir:pathlib.Path, label:str) -> pathlib.Path:
	'''
	Return the path of an upstream source tree, downloading it unless one was supplied.

	:param tag: the upstream git tag
	:param override: a path given on the command line, or None
	:param workdir: a scratch directory to download into
	:param label: "target" or "previous", for messages
	:returns: the path of the source tree
	'''
	if override is not None:
		path = pathlib.Path(override).resolve()
		if not path.is_dir():
			raise SyncError(f"The {label} upstream directory '{path}' does not exist.")
		print(f"  using the {label} upstream tree at {path}")
		return path
	return download_upstream(tag, workdir / label)


def read_upstream_file(tree:pathlib.Path, filename:str) -> str:
	'''
	Read one file from an unpacked upstream source tree.

	:param tree: the unpacked upstream source tree
	:param filename: the base name of the file to read
	:returns: the file's contents
	'''
	path = tree / filename
	if not path.is_file():
		raise SyncError(f"The upstream tree at '{tree}' has no file '{filename}'. Either the release "
						f"reorganized its sources or the wrong directory was given; the vendored file list "
						f"in this script needs a human's attention.")
	return path.read_text()


# ------------------------------------------------------------------
# Running the gates
# ------------------------------------------------------------------

def run_command(command:list, cwd:pathlib.Path):
	'''
	Run a subprocess, raising a :class:`SyncError` if it fails.

	:param command: the command and its arguments
	:param cwd: the working directory to run it in
	'''
	print(f"  $ {' '.join(command)}   (in {cwd})")
	result = subprocess.run(command, cwd=str(cwd))
	if result.returncode != 0:
		raise SyncError(f"Command failed with exit status {result.returncode}: {' '.join(command)}\n"
						f"The vendored sources have already been written; inspect the failure, then either fix "
						f"it or 'git checkout' the C sources and setup.py to abandon the sync.")


def run_gates():
	'''
	Reinstall the package and run the guard tests followed by the full suite.

	The tests import the installed package, not the working tree, so a sync that is
	not reinstalled first would test the old C code and pass for the wrong reason.
	'''
	print("\nReinstalling the package (the tests import the installed package, not the working tree):")
	run_command([sys.executable, "-m", "pip", "install", "--force-reinstall", "--no-deps", "."], SOURCE_DIR)
	print("\nRunning the source invariant guard tests:")
	run_command([sys.executable, "-m", "pytest", "-q", "test_source_invariants.py"], TESTS_DIR)
	print("\nRunning the full test suite:")
	run_command([sys.executable, "-m", "pytest", "-q"], TESTS_DIR)


# ------------------------------------------------------------------
# Orchestration
# ------------------------------------------------------------------

def check_no_undocumented_divergence(previous_tree:pathlib.Path):
	'''
	Confirm the vendored files are the previous upstream release plus exactly the
	patches this script knows about.

	Anything else in them is a local change nobody recorded, and the sync would
	overwrite it without a trace.

	:param previous_tree: the unpacked upstream tree the vendored files came from
	'''
	for filename in VENDORED_FILES:
		vendored = (VENDOR_DIR / filename).read_text()
		for patch in patches_for(filename):
			vendored = unapply_patch(patch, vendored)
		upstream = read_upstream_file(previous_tree, filename)
		if vendored != upstream:
			diff = "\n".join(list(difflib.unified_diff(upstream.splitlines(), vendored.splitlines(),
													   fromfile=f"upstream/{filename}",
													   tofile=f"vendored/{filename} (known patches reversed)",
													   lineterm=""))[:60])
			raise SyncError(f"The vendored '{filename}' is not the previous upstream release plus the patches "
							f"this script knows about. There is an undocumented local divergence, and syncing "
							f"would silently discard it.\n"
							f"WHAT TO DO: identify the change below, then either add it to PATCHES in this "
							f"script (and to CLAUDE.md) or revert it.\n\n{diff}")
	print("  the vendored files are the previous upstream release plus the known patches: OK")


def summarize_change(filename:str, old_text:str, new_text:str) -> str:
	'''
	Describe how a vendored file would change, in one line.

	:param filename: the file's base name
	:param old_text: the current contents
	:param new_text: the contents that would be written
	:returns: a one-line description
	'''
	if old_text == new_text:
		return f"  {filename}: unchanged"
	added = removed = 0
	for line in difflib.unified_diff(old_text.splitlines(), new_text.splitlines(), lineterm="", n=0):
		if line.startswith("+") and not line.startswith("+++"):
			added += 1
		elif line.startswith("-") and not line.startswith("---"):
			removed += 1
	return f"  {filename}: {added} line(s) added, {removed} removed"


def sync(tag:str, from_tag:str=None, dry_run:bool=False, skip_tests:bool=False,
		 upstream_dir:str=None, previous_dir:str=None) -> int:
	'''
	Perform the sync. See the module docstring for the sequence.

	:param tag: the upstream tag to sync to
	:param from_tag: the upstream tag the vendored files currently come from; derived from setup.py if None
	:param dry_run: report what would change and write nothing
	:param skip_tests: do not reinstall the package and run pytest
	:param upstream_dir: use this unpacked upstream tree instead of downloading `tag`
	:param previous_dir: use this unpacked upstream tree instead of downloading `from_tag`
	:returns: a process exit status
	'''
	setup_text = SETUP_PATH.read_text()
	current_version = setup_q3c_version(setup_text)
	new_version = version_from_tag(tag)
	if from_tag is None:
		from_tag = f"v{current_version}"

	print(f"Syncing the vendored Q3C sources from {from_tag} to {tag}"
		  + (" (DRY RUN - nothing will be written)" if dry_run else ""))

	with tempfile.TemporaryDirectory(prefix="qlsc-q3c-sync-") as tmp:
		workdir = pathlib.Path(tmp)
		(workdir / "target").mkdir()
		(workdir / "previous").mkdir()

		print("\n1. Fetching upstream:")
		target_tree = upstream_tree(tag, upstream_dir, workdir, "target")
		if previous_dir is None and from_tag == tag:
			# re-syncing the current release against itself, the end-to-end self check
			print(f"  the previous release is the same tag; reusing the tree just fetched")
			previous_tree = target_tree
		else:
			previous_tree = upstream_tree(from_tag, previous_dir, workdir, "previous")

		print("\n2. Checking the vendored files against the release they came from:")
		check_no_undocumented_divergence(previous_tree)

		print("\n3. Deciding on the local patches:")
		new_contents = {}
		patch_states = {}
		for filename in VENDORED_FILES:
			text = read_upstream_file(target_tree, filename)
			for patch in patches_for(filename):
				text, state = apply_patch(patch, text)
				patch_states[patch.name] = state
				if state is PatchState.NEEDED:
					print(f"  {patch.name}: still needed upstream, re-applied "
						  f"({patch.occurrences} site(s))")
				else:
					print(f"  {patch.name}: PRESENT UPSTREAM, patch dropped")
			new_contents[filename] = text

		print("\n4. Checking the C interface the wrapper depends on:")
		differences = interface_differences(read_upstream_file(previous_tree, HEADER_FILE),
											read_upstream_file(target_tree, HEADER_FILE),
											WRAPPER_PATH.read_text())
		if differences:
			raise SyncError("The C interface used by qlsc_c_module.c changed between "
							f"{from_tag} and {tag}:\n\n" + "\n".join(differences)
							+ "\n\nWHAT TO DO: update qlsc_c_module.c for the new interface by hand, then "
							  "re-run this script.")
		print(f"  every {HEADER_FILE} declaration, struct and macro the wrapper uses is unchanged: OK")

		print("\n5. Version string:")
		upstream_control = control_version(read_upstream_file(target_tree, "q3c.control"))
		if upstream_control != new_version:
			raise SyncError(f"The tag '{tag}' implies Q3C_VERSION '{new_version}', but upstream's q3c.control "
							f"at that tag says '{upstream_control}'. Resolve the discrepancy by hand.")
		print(f"  setup.py Q3C_VERSION: {current_version} -> {new_version} (confirmed by upstream's q3c.control)")

		print("\n6. Changes to the vendored files:")
		changed = False
		for filename in VENDORED_FILES:
			old_text = (VENDOR_DIR / filename).read_text()
			print(summarize_change(filename, old_text, new_contents[filename]))
			changed = changed or old_text != new_contents[filename]

		if dry_run:
			print("\nDry run: no files were written, nothing was installed and no tests were run.")
			_print_checklist(tag, patch_states, changed or current_version != new_version)
			return 0

		print("\n7. Writing:")
		for filename in VENDORED_FILES:
			(VENDOR_DIR / filename).write_text(new_contents[filename])
			print(f"  wrote {VENDOR_DIR / filename}")
		if current_version != new_version:
			SETUP_PATH.write_text(set_setup_q3c_version(setup_text, new_version))
			print(f"  wrote {SETUP_PATH} (Q3C_VERSION = {new_version})")

	if skip_tests:
		print("\n--skip-tests given: the package was NOT reinstalled and no tests were run. "
			  "Do both before committing.")
	else:
		run_gates()

	_print_checklist(tag, patch_states, True)
	return 0


def _print_checklist(tag:str, patch_states:dict, changed:bool):
	'''
	Print the human work this script cannot do.

	:param tag: the upstream tag synced to
	:param patch_states: each patch name mapped to its :class:`PatchState`
	:param changed: whether anything actually changed
	'''
	print("\n" + "=" * 78)
	if not changed:
		print("Nothing changed: the vendored sources already match this release.")
		print("=" * 78)
		return
	print("STILL TO DO BY HAND:")
	print(f"  * Re-read the open issues at {UPSTREAM_ISSUES_URL} and check whether any bug that was")
	print("    dormant in qlsc (unreachable from the wrapper) is now reachable, and whether any that")
	print("    were reachable have been fixed. CLAUDE.md lists what was assessed last time.")
	print("  * Check whether the release added a source file the extension now needs; this script")
	print(f"    only ever copies {', '.join(VENDORED_FILES)}.")
	for name, state in patch_states.items():
		if state is PatchState.ALREADY_UPSTREAM:
			patch = next(p for p in PATCHES if p.name == name)
			print(f"  * The patch '{name}' is now upstream. Remove it from PATCHES in this script, remove")
			print(f"    the guard test '{patch.guard_test}' from source/tests/test_source_invariants.py")
			print("    (which will FAIL until you do), and update the CLAUDE.md sync section.")
	print(f"  * Update CLAUDE.md (the vendored Q3C version is stated there) and CHANGELOG.md for {tag}.")
	print("=" * 78)


def main(argv=None) -> int:
	'''
	Command line entry point.

	:param argv: argument list, defaulting to sys.argv[1:]
	:returns: a process exit status
	'''
	parser = argparse.ArgumentParser(description="Sync the vendored Q3C C sources to an upstream release.",
									 epilog="See the docstring at the top of this file for details.")
	parser.add_argument("tag", help="the upstream q3c tag to sync to, e.g. v2.0.6")
	parser.add_argument("--from-tag", default=None,
						help="the upstream tag the vendored files currently come from "
							 "(default: derived from the Q3C_VERSION in setup.py)")
	parser.add_argument("--dry-run", action="store_true",
						help="report what would change; write nothing, install nothing, run no tests")
	parser.add_argument("--skip-tests", action="store_true",
						help="do not reinstall the package and run pytest afterwards")
	parser.add_argument("--upstream-dir", default=None,
						help="use this already-unpacked upstream tree instead of downloading the tag")
	parser.add_argument("--previous-dir", default=None,
						help="use this already-unpacked upstream tree instead of downloading --from-tag")
	args = parser.parse_args(argv)

	try:
		return sync(tag=args.tag,
					from_tag=args.from_tag,
					dry_run=args.dry_run,
					skip_tests=args.skip_tests,
					upstream_dir=args.upstream_dir,
					previous_dir=args.previous_dir)
	except SyncError as e:
		print(f"\nSYNC STOPPED\n\n{e}\n", file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())
