
import ast
import pathlib
import re

import pytest

# A gate over the whole tree: nothing shipped or run may use syntax or standard
# library modules newer than the Python version setup.py declares support for.
#
# This exists because the invariant is mechanically decidable and was broken once:
# scripts/sync_q3c.py used dataclasses (3.7) while python_requires said >=3.6, and
# because the test suite imports that script the whole suite would have failed on a
# supported interpreter. Fixing that one file would have left the rule unenforced,
# so the rule is enforced here instead.
#
# Blind spots, stated rather than assumed covered:
#  * only syntax and stdlib module names are checked, not new METHODS or keyword
#    arguments on existing objects (subprocess.run(text=...), 3.7, or
#    zip(strict=...), 3.10, would pass this gate);
#  * f-string self-documenting expressions, f"{x=}" (3.8), are not detected - they
#    do not survive into a distinguishable AST node across versions;
#  * third-party APIs (NumPy, astropy) are not version checked here.
# The reliable check for those remains running the suite on the oldest supported
# interpreter.

REPO_SOURCE = pathlib.Path(__file__).resolve().parent.parent
REPO_ROOT = REPO_SOURCE.parent
SETUP_PATH = REPO_SOURCE / "setup.py"

# module name -> the Python version that first shipped it in the standard library
STDLIB_MODULE_VERSIONS = {
	"dataclasses": (3, 7),
	"contextvars": (3, 7),
	"importlib.resources": (3, 7),
	"graphlib": (3, 9),
	"zoneinfo": (3, 9),
	"tomllib": (3, 11)
}

def python_floor() -> tuple:
	'''
	Return the minimum Python version setup.py declares support for, as a tuple.
	'''
	match = re.search(r"python_requires\s*=\s*['\"]>=\s*(\d+)\.(\d+)['\"]", SETUP_PATH.read_text(encoding="utf-8"))
	assert match is not None, "could not read python_requires from setup.py"
	return (int(match.group(1)), int(match.group(2)))

def python_files() -> list:
	'''
	Return every Python file in the repository that is shipped or run: the package,
	the repository scripts, and the tests themselves.
	'''
	files = sorted((REPO_SOURCE / "qlsc").rglob("*.py"))
	files.append(SETUP_PATH)
	files += sorted((REPO_SOURCE / "tests").glob("*.py"))
	scripts = REPO_ROOT / "scripts"
	if scripts.is_dir():
		files += sorted(scripts.glob("*.py"))
	return files

def annotation_nodes(tree:ast.AST):
	'''
	Yield every AST node that appears inside a type annotation.
	'''
	for node in ast.walk(tree):
		annotations = []
		if isinstance(node, ast.AnnAssign) and node.annotation is not None:
			annotations.append(node.annotation)
		elif isinstance(node, ast.arg) and node.annotation is not None:
			annotations.append(node.annotation)
		elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.returns is not None:
			annotations.append(node.returns)
		for annotation in annotations:
			for inner in ast.walk(annotation):
				yield inner

def version_violations(path:pathlib.Path, floor:tuple) -> list:
	'''
	Return a list of (line, required version, description) for every construct in the
	file that needs a Python newer than `floor`.
	'''
	tree = ast.parse(path.read_text(encoding="utf-8"), str(path))
	found = []

	def report(node, version, description):
		if version > floor:
			found.append((getattr(node, "lineno", 0), version, description))

	for node in ast.walk(tree):
		if isinstance(node, ast.NamedExpr):
			report(node, (3, 8), "assignment expression ':='")
		elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and getattr(node.args, "posonlyargs", []):
			report(node, (3, 8), "positional-only parameters")
		elif node.__class__.__name__ == "Match":
			report(node, (3, 10), "match statement")
		elif isinstance(node, ast.ImportFrom) and node.module == "__future__":
			for alias in node.names:
				if alias.name == "annotations":
					report(node, (3, 7), "from __future__ import annotations")
		elif isinstance(node, ast.Import):
			for alias in node.names:
				if alias.name in STDLIB_MODULE_VERSIONS:
					report(node, STDLIB_MODULE_VERSIONS[alias.name], f"import {alias.name}")
		elif isinstance(node, ast.ImportFrom) and node.module in STDLIB_MODULE_VERSIONS:
			report(node, STDLIB_MODULE_VERSIONS[node.module], f"from {node.module} import ...")

	# subscripted builtin generics (list[int]) are 3.9; only annotations are inspected,
	# so an ordinary subscript of a variable that happens to be named "list" is not flagged
	for node in annotation_nodes(tree):
		if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name) and \
		   node.value.id in ("list", "dict", "set", "tuple", "frozenset", "type"):
			report(node, (3, 9), f"subscripted builtin '{node.value.id}[...]' in an annotation")

	return found

def test_python_requires_is_readable():
	'''
	Test that the declared floor can be read; the rest of this file depends on it.
	'''
	floor = python_floor()
	assert floor >= (3, 0)

def test_no_construct_newer_than_the_declared_floor():
	'''
	Test that no shipped or executed Python file uses syntax or a standard library
	module newer than setup.py's python_requires.
	'''
	floor = python_floor()
	problems = []
	for path in python_files():
		for lineno, version, description in version_violations(path, floor):
			problems.append(f"{path.relative_to(REPO_ROOT)}:{lineno}: {description} needs Python "
							f"{version[0]}.{version[1]}, but python_requires is >={floor[0]}.{floor[1]}")
	assert not problems, ("Python constructs newer than the declared floor:\n  " + "\n  ".join(problems)
						  + "\nEither rewrite them or raise python_requires in setup.py deliberately.")

def test_the_gate_detects_a_violation():
	'''
	Test that the check above actually fires; a gate that cannot fail guards nothing.
	'''
	sample = pathlib.Path(__file__).parent / "__version_floor_sample__.py"
	sample.write_text("import dataclasses\n\ndef f(values: list[int]) -> dict[str,int]:\n\treturn {}\n", encoding="utf-8")
	try:
		violations = version_violations(sample, (3, 6))
		descriptions = [d for _, _, d in violations]
		assert "import dataclasses" in descriptions
		assert any("subscripted builtin" in d for d in descriptions)
		# and nothing is reported when the floor is high enough
		assert version_violations(sample, (3, 12)) == []
	finally:
		sample.unlink()


# A second gate over the same files: text I/O must name its encoding. Without one,
# Python uses the locale's encoding, which is cp1252 on Windows, and a file holding
# any non-ASCII character (the README, CLAUDE.md) fails to decode there. This broke
# the Windows wheel tests once (test_readme_appendix.py), with some thirty more
# instances in the tree waiting behind it, so the rule is enforced here instead.
#
# Blind spots, stated rather than assumed covered:
#  * pathlib.Path.open(), io.open() and codecs.open() are not checked (an attribute
#    named "open" is too often tarfile.open or gzip.open to judge by name);
#  * subprocess's text=True decodes child output in the locale's encoding;
#  * an open() whose mode is not a literal string is reported, conservatively;
#  * the file set is python_files(): docs/*.py (check_example_outputs.py, conf.py)
#    and .github/scripts/ are not gated - they run only on Linux CI.

def unencoded_text_io(path:pathlib.Path) -> list:
	'''
	Return a list of (line, description) for every text-mode read or write in the
	file that does not pass an explicit encoding.
	'''
	tree = ast.parse(path.read_text(encoding="utf-8"), str(path))
	found = []
	for node in ast.walk(tree):
		if not isinstance(node, ast.Call):
			continue
		if any(keyword.arg == "encoding" for keyword in node.keywords):
			continue
		function = node.func
		if isinstance(function, ast.Attribute) and function.attr in ("read_text", "write_text"):
			# read_text(encoding) / write_text(data, encoding) may pass it positionally
			positional_encoding = len(node.args) >= (1 if function.attr == "read_text" else 2)
			if not positional_encoding:
				found.append((node.lineno, f".{function.attr}() without an encoding"))
		elif isinstance(function, ast.Name) and function.id == "open":
			mode = node.args[1] if len(node.args) >= 2 else next(
				(keyword.value for keyword in node.keywords if keyword.arg == "mode"), None)
			if isinstance(mode, ast.Constant) and isinstance(mode.value, str) and "b" in mode.value:
				continue # binary mode has no encoding
			if len(node.args) >= 4:
				continue # encoding given positionally: open(file, mode, buffering, encoding)
			found.append((node.lineno, "open() in text mode without an encoding"))
	return found

def test_text_io_names_its_encoding():
	'''
	Test that no shipped or run Python file reads or writes text in the locale's
	default encoding.
	'''
	offences = []
	for path in python_files():
		for line, description in unencoded_text_io(path):
			offences.append(f"{path.relative_to(REPO_ROOT)}:{line}: {description}")
	assert not offences, ("text I/O without an explicit encoding (cp1252 on Windows):\n"
						  + "\n".join(offences) + "\nPass encoding=\"utf-8\".")
