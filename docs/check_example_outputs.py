#!/usr/bin/env python
'''
Re-run the code blocks in the documentation and check that they still print what
the pages say they print.

The documentation makes a promise — that every output block was produced by
running the code above it — and a promise that nothing enforces decays. This
script enforces it.

WHAT IT CHECKS
--------------
A **verified pair** is a ``.. code-block:: python`` immediately followed by a
``.. code-block:: text``. Those pairs, and only those, are executed: the Python
blocks of one page run in order in a single namespace (pages build state up
across blocks), each one's standard output is captured, and it is compared
against the text block that follows it. Each page runs in its own empty
temporary directory, so blocks that write files are hermetic.

Output is compared after normalizing trailing whitespace, leading and trailing blank lines,
hexadecimal object addresses (``0x10619ad50`` in a ``repr``), and the block's
common leading indentation — docutils strips the last of those when it renders a
literal block, so reST cannot represent output whose every line starts with the
same spaces. Indentation *relative* to the block is still compared, as is
everything else.

SKIP MARKERS
------------
A block that cannot run hermetically is marked in the page source itself, with a
reST comment on the line before its ``code-block`` directive::

	.. check-outputs: run-only -- wall-clock timings differ between machines

	.. check-outputs: skip -- needs network access to VizieR

``run-only`` executes the block (so later blocks still have its state) but does
not compare its output. ``skip`` does neither. ``setup`` is for a block that has
no output block of its own but whose imports or objects the blocks below it
use — a tutorial page that introduces its imports in their own block, for
instance; it is executed and nothing is compared. A whole page carries::

	.. check-outputs-page: skip -- reason it can never run here

	.. check-outputs-page: needs-network -- fetches the Bright Star Catalogue from VizieR

``skip`` never runs the page. ``needs-network`` runs it only when ``--network`` is
given, so the default run stays hermetic and offline while the page can still be
verified deliberately; ``needs-cornish`` does the same for ``--cornish`` and the
optional ``cornish`` package, which is not installed with qlsc. Both are listed in
``OPTIONAL_PAGE_MODES``; a page carrying one is *reported* as skipped, with its
reason, rather than passing silently.

Every marker must carry a reason after ``--``; the script refuses a marker
without one. The reasons live in the page source, where a reader editing the page
sees them, rather than in a list kept somewhere else.

BLIND SPOTS — what this does NOT check
--------------------------------------
* **Unpaired Python blocks are not run** unless marked ``setup``. A
  ``code-block:: python`` with no text block after it and no marker is invisible to
  this script: it is never executed, so it can rot into something that would not
  even import. Many are deliberately illustrative fragments with undefined names,
  so running them is not possible; the summary reports how many there are per page
  so the number cannot grow silently.
* **Prose is not checked.** A number quoted in a sentence rather than in an output
  block ("about 51,000 per square degree") is not compared against anything. Nor
  are the figures, whose numbers come from ``generate_figures.py`` and must agree
  with the pages by hand.
* **Skipped and run-only blocks' outputs are never verified**, by construction.
  Their pasted output is a standing manual duty; ``--list`` prints them.
* **Only ``.. code-block::`` is recognized**, not ``::`` literal blocks or
  ``doctest`` blocks. The documentation uses ``code-block`` throughout; a page
  that switched style would be silently uncovered, so the summary's total pair
  count is worth watching.
* Text blocks that are *not* program output (SQL, shell transcripts, formulas)
  are only examined when they directly follow a Python block, so a text block
  after a ``code-block:: sql`` or after prose is correctly ignored.

USAGE
-----
::

	cd docs
	python check_example_outputs.py            # every page
	python check_example_outputs.py --list     # what is covered and what is skipped
	python check_example_outputs.py examples/cone_search.rst

Exits non-zero if any pair mismatches or raises. Needs the same packages the
examples do: numpy, astropy and qlsc itself.
'''

import argparse
import contextlib
import io
import os
import re
import sys
import tempfile
import traceback

DOCS_ROOT = os.path.dirname(os.path.abspath(__file__))

BLOCK_RE = re.compile(r'^(\s*)\.\.\s+code-block::\s*(\S+)\s*$')
MARKER_RE = re.compile(r'^\s*\.\.\s+check-outputs:\s*(\S+)\s*(?:--\s*(.*))?$')
PAGE_MARKER_RE = re.compile(r'^\s*\.\.\s+check-outputs-page:\s*(\S+)\s*(?:--\s*(.*))?$')
ADDRESS_RE = re.compile(r'0x[0-9a-fA-F]+')

MODES = ("skip", "run-only", "setup")

#: Page modes that run the page only when the matching command line flag is given,
#: mapped to that flag and its help text. A page whose every block needs something the
#: default run does not have - a network service, an optional package - carries one of
#: these, and is reported as skipped with its reason rather than passing quietly.
#:
#: This table is the only place a mode is declared: main() builds the command line
#: options from it and reads the results back out of it, so adding a third mode is one
#: entry here and a marker in the page, never a second edit somewhere that has to agree
#: with the first.
OPTIONAL_PAGE_MODES = {
	"needs-network": ("--network", "also run pages marked 'needs-network' (they fetch from "
								   "a remote service)"),
	"needs-cornish": ("--cornish", "also run pages marked 'needs-cornish' (they need the "
								   "optional 'cornish' package, which is not installed with qlsc)"),
	}
PAGE_MODES = ("skip",) + tuple(OPTIONAL_PAGE_MODES)

def _option_destination(flag:str) -> str:
	'''
	Return the attribute argparse stores a long option under, derived the way argparse
	derives it: the leading dashes are dropped and any remaining dash becomes an
	underscore. Written out rather than assumed so that a mode whose flag has a dash
	inside it cannot silently look up an attribute that is not there.

	:param flag: the long option, e.g. "--needs-network"
	'''
	return flag.lstrip("-").replace("-", "_")


class Block:
	def __init__(self, language, body, line, indent, body_start, body_end, mode=None, reason=None):
		self.language = language
		self.body = body
		self.line = line
		self.indent = indent          # the directive's own indentation
		self.body_start = body_start  # 0-based index of the first body line
		self.body_end = body_end      # 0-based index one past the last body line
		self.mode = mode
		self.reason = reason
		self.expected = None          # the text block that follows, if any

	@property
	def paired(self):
		return self.expected is not None


def _dedent(lines):
	'''Strip the common leading whitespace from a block body, tabs included.'''
	widths = [len(l) - len(l.lstrip()) for l in lines if l.strip()]
	if not widths:
		return []
	cut = min(widths)
	return [l[cut:] if l.strip() else "" for l in lines]


def parse_page(path):
	'''
	Return (blocks, page_skip_reason) for one reST file.

	Blocks are returned in source order with their pairing already resolved.
	'''
	with open(path, encoding="utf-8") as f:
		lines = f.read().splitlines()

	page_mode = page_skip = None
	for line in lines:
		m = PAGE_MARKER_RE.match(line)
		if m:
			mode, reason = m.group(1), (m.group(2) or "").strip()
			if mode not in PAGE_MODES:
				raise ValueError(f"{path}: check-outputs-page understands "
								 f"{' and '.join(PAGE_MODES)}, got '{mode}'")
			if not reason:
				raise ValueError(f"{path}: the check-outputs-page marker must give a reason after '--'")
			page_mode, page_skip = mode, reason
			break

	blocks = []
	pending = None            # (mode, reason) from a marker seen just above
	i = 0
	while i < len(lines):
		marker = MARKER_RE.match(lines[i])
		if marker:
			mode, reason = marker.group(1), (marker.group(2) or "").strip()
			if mode not in MODES:
				raise ValueError(f"{path}:{i+1}: unknown check-outputs mode '{mode}'; "
								 f"expected one of {', '.join(MODES)}")
			if not reason:
				raise ValueError(f"{path}:{i+1}: the check-outputs marker must give a reason after '--'")
			pending = (mode, reason)
			i += 1
			continue

		m = BLOCK_RE.match(lines[i])
		if not m:
			# a marker only applies to the directive that follows it immediately,
			# ignoring blank lines; anything else in between drops it
			if pending and lines[i].strip():
				pending = None
			i += 1
			continue

		indent, language = m.group(1), m.group(2)
		start_line = i + 1
		i += 1
		# skip the directive's options and the blank line before its body
		while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith(":")):
			i += 1
		body = []
		body_start = i
		while i < len(lines):
			line = lines[i]
			if not line.strip():
				body.append("")
				i += 1
				continue
			if len(line) - len(line.lstrip()) <= len(indent):
				break
			body.append(line)
			i += 1

		mode, reason = pending if pending else (None, None)
		pending = None
		blocks.append(Block(language, "\n".join(_dedent(body)).strip("\n"),
							start_line, indent, body_start, i, mode, reason))

	# pair each python block with the text block that immediately follows it
	for a, b in zip(blocks, blocks[1:]):
		if a.language == "python" and b.language == "text":
			a.expected = b

	return [b for b in blocks if b.language == "python"], (page_mode, page_skip)


def normalize(text):
	'''Remove the differences that are not the page's fault.'''
	text = ADDRESS_RE.sub("0xADDRESS", text)
	lines = [l.rstrip() for l in text.splitlines()]
	while lines and not lines[-1]:
		lines.pop()
	while lines and not lines[0]:
		lines.pop(0)
	# docutils strips a literal block's common indentation, so a reST text block
	# cannot represent output whose every line begins with the same spaces. Both
	# sides are dedented so the comparison tests what the reader actually sees;
	# indentation *relative* to the block is still compared.
	return "\n".join(_dedent(lines))


def check_page(path, verbose=False, enabled_modes=()):
	'''
	Run one page. Returns (failures, checked, run_only, skipped, unpaired).

	:param path: the reST file to run
	:param verbose: report every block, not only the failures
	:param enabled_modes: the :py:data:`OPTIONAL_PAGE_MODES` the caller asked for; a page carrying one that is not here is skipped with its reason
	'''
	relative = os.path.relpath(path, DOCS_ROOT)
	blocks, (page_mode, page_skip) = parse_page(path)

	if page_mode == "skip" or (page_mode in OPTIONAL_PAGE_MODES and page_mode not in enabled_modes):
		suffix = f" (pass {OPTIONAL_PAGE_MODES[page_mode][0]} to run it)" if page_mode in OPTIONAL_PAGE_MODES else ""
		print(f"SKIP  {relative}: {page_skip}{suffix}")
		return 0, 0, 0, len(blocks), 0

	paired = [b for b in blocks if b.paired]
	runnable = [b for b in blocks if b.paired or b.mode == "setup"]
	if not paired:
		return 0, 0, 0, 0, len(blocks)

	failures = checked = run_only = skipped = 0
	unpaired = len([b for b in blocks if not b.paired and b.mode != "setup"])
	namespace = {"__name__": "__docs__"}
	previous_directory = os.getcwd()

	with tempfile.TemporaryDirectory(prefix="qlsc-docs-check-") as working_directory:
		os.chdir(working_directory)
		try:
			for block in runnable:
				where = f"{relative}:{block.line}"
				if block.mode == "skip":
					skipped += 1
					if verbose:
						print(f"  skip     {where}: {block.reason}")
					continue

				captured = io.StringIO()
				try:
					with contextlib.redirect_stdout(captured):
						exec(compile(block.body, where, "exec"), namespace)
				except Exception:
					failures += 1
					print(f"FAIL  {where}: the block raised")
					print(textwrap_indent(traceback.format_exc().rstrip(), "        "))
					continue

				if block.mode in ("run-only", "setup"):
					run_only += 1
					if verbose:
						print(f"  run-only {where}: {block.reason}")
					continue

				got = normalize(captured.getvalue())
				want = normalize(block.expected.body)
				checked += 1
				if got != want:
					failures += 1
					print(f"FAIL  {where}: output does not match the page")
					print(textwrap_indent("--- page says ---\n" + want, "        "))
					print(textwrap_indent("--- got ---\n" + got, "        "))
				elif verbose:
					print(f"  ok       {where}")
		finally:
			os.chdir(previous_directory)

	status = "FAIL" if failures else "ok  "
	print(f"{status}  {relative}: {checked} checked, {run_only} run-only, "
		  f"{skipped} skipped, {unpaired} unpaired")
	return failures, checked, run_only, skipped, unpaired


def textwrap_indent(text, prefix):
	return "\n".join(prefix + l if l else l for l in text.splitlines())


def rst_files(names):
	if names:
		return [os.path.join(DOCS_ROOT, n) for n in names]
	found = []
	for root, directories, files in os.walk(DOCS_ROOT):
		directories[:] = [d for d in directories if d not in ("_build", "_static", "_templates")]
		found.extend(os.path.join(root, f) for f in files if f.endswith(".rst"))
	return sorted(found)


def list_coverage(paths):
	for path in paths:
		relative = os.path.relpath(path, DOCS_ROOT)
		blocks, (page_mode, page_skip) = parse_page(path)
		paired = [b for b in blocks if b.paired]
		setups = [b for b in blocks if not b.paired and b.mode == "setup"]
		if page_mode:
			print(f"{relative}: PAGE {page_mode.upper()} -- {page_skip}")
			continue
		if not blocks:
			continue
		print(f"{relative}: {len(paired)} verified pairs, "
			  f"{len(blocks)-len(paired)-len(setups)} unpaired blocks")
		for block in blocks:
			if block.mode:
				print(f"    line {block.line}: {block.mode} -- {block.reason}")


def main():
	parser = argparse.ArgumentParser(description="Re-run the documentation's example code and "
												 "diff it against the output the pages show.")
	parser.add_argument("pages", nargs="*", help="reST files to check, relative to docs/; "
												"the default is all of them")
	parser.add_argument("--list", action="store_true", help="report coverage and skips, run nothing")
	parser.add_argument("-v", "--verbose", action="store_true", help="report every block")
	for flag, help_text in OPTIONAL_PAGE_MODES.values():
		parser.add_argument(flag, action="store_true", help=help_text)
	args = parser.parse_args()

	enabled_modes = tuple(mode for mode, (flag, _) in OPTIONAL_PAGE_MODES.items()
						  if getattr(args, _option_destination(flag)))

	paths = rst_files(args.pages)
	if args.list:
		list_coverage(paths)
		return 0

	totals = [0, 0, 0, 0, 0]
	for path in paths:
		result = check_page(path, verbose=args.verbose, enabled_modes=enabled_modes)
		totals = [t + r for t, r in zip(totals, result)]

	failures, checked, run_only, skipped, unpaired = totals
	print()
	print(f"{checked} blocks checked, {run_only} run-only, {skipped} skipped, "
		  f"{unpaired} unpaired (never executed), {failures} failures")
	if failures:
		print("Some documented output no longer matches. Re-run the code and paste what it "
			  "really prints; never edit the expected output by hand.")
	return 1 if failures else 0


if __name__ == "__main__":
	sys.exit(main())
