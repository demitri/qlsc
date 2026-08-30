
import pathlib
import re

# Mechanical guards over the C source itself, for invariants that the runtime
# tests cannot see failing (a reintroduced stack overflow is a silent 32-byte
# scribble, not a crash).

C_CODE_DIR = pathlib.Path(__file__).parent.parent / "qlsc" / "c_code"

def test_q3c_poly_points_buffer_patch():
	'''
	Test that the local points[8] patch to the vendored q3c_poly.c is present.

	Upstream q3c v2.0.5 declares points[4] in q3c_check_sphere_point_in_poly, an
	ASan-confirmed stack buffer overflow reachable through QLSC.point_in_polygon
	(q3c_multi_face_check writes up to four coordinate pairs). A future upstream
	sync that copies the file verbatim would silently reintroduce it; this test
	fails loudly instead. See CLAUDE.md. Remove this test only when upstream
	ships the fix and the sync drops the local patch.
	'''
	source = (C_CODE_DIR / "q3c" / "q3c_poly.c").read_text()
	assert "q3c_coord_t points[8];" in source, "the local points[8] patch to q3c_poly.c has been lost (probably by an upstream sync)"
	assert "q3c_coord_t points[4];" not in source, "upstream's undersized points[4] buffer has been reintroduced into q3c_poly.c"

def test_q3cube_res_depth_patch():
	'''
	Test that the local subdivision-depth patch to the vendored q3cube.c is present.

	Upstream q3c v2.0.5 computes the quadtree refinement depth as the ratio
	nside/n0 instead of its log2 (upstream issue #51), overflowing a signed int
	(2 * 2^30) for tiny query regions - UBSan-confirmed, producing garbage
	coordinates. Remove this test only when upstream ships the fix and the sync
	drops the local patch. See CLAUDE.md.
	'''
	source = (C_CODE_DIR / "q3c" / "q3cube.c").read_text()
	assert "res_depth = nside / n0;" not in source, "upstream's overflowing res_depth ratio has been reintroduced into q3cube.c"
	assert source.count("subdiv_ratio") >= 3, "the local log2 subdivision-depth patch to q3cube.c has been lost (probably by an upstream sync)"

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
