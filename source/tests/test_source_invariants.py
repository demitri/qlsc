
import pathlib
import re

# Mechanical guards over the C source itself, for invariants that the runtime
# tests cannot see failing.
#
# The vendored q3c files carry no local patches (q3c v2.0.6 fixed upstream the
# two bugs qlsc used to patch, #50 and #51), so there is no patch guard here:
# scripts/sync_q3c.py's divergence check refuses any vendored difference not
# registered in its PATCHES, and the runtime regressions for both bugs live in
# test_point_in_polygon.py and test_radial_query_center_coverage.py.

C_CODE_DIR = pathlib.Path(__file__).parent.parent / "qlsc" / "c_code"

def test_capsule_pointer_null_checks():
	'''
	Test that every PyCapsule_GetPointer call in the wrapper is followed by a
	NULL check; an unchecked one segfaults the interpreter when handed an
	invalid capsule.
	'''
	# Note: this is a whole-file count, not a per-site check - it cannot say WHICH
	# site is unguarded, and guards clustered elsewhere would satisfy it. It is a
	# smoke guard; the per-function runtime check lives in test_invalid_inputs.py.
	source = (C_CODE_DIR / "qlsc_c_module.c").read_text()
	# strip line comments so commented-out code is not counted
	source = re.sub(r"//[^\n]*", "", source)
	calls = source.count("PyCapsule_GetPointer")
	guards = source.count("if (hprm == NULL)")
	assert calls > 0
	assert guards == calls, f"{calls} PyCapsule_GetPointer call(s) but {guards} NULL guard(s) in qlsc_c_module.c; every call must be immediately followed by a NULL check"
