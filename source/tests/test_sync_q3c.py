
import importlib.util
import pathlib

import pytest

# Tests for the upstream sync script, scripts/sync_q3c.py.
#
# The script lives in the repository, not in the distributed package, so it is not
# present when the tests are run from an unpacked sdist. That case is skipped
# explicitly (and says so) rather than being quietly passed over.

SCRIPT_PATH = pathlib.Path(__file__).resolve().parent.parent.parent / "scripts" / "sync_q3c.py"

if not SCRIPT_PATH.is_file():
	pytest.skip(f"scripts/sync_q3c.py is not present at {SCRIPT_PATH}; it ships in the repository, "
				"not in the qlsc package.", allow_module_level=True)

_spec = importlib.util.spec_from_file_location("qlsc_sync_q3c", SCRIPT_PATH)
sync_q3c = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sync_q3c)

VENDOR_DIR = pathlib.Path(__file__).resolve().parent.parent / "qlsc" / "c_code" / "q3c"
WRAPPER_PATH = pathlib.Path(__file__).resolve().parent.parent / "qlsc" / "c_code" / "qlsc_c_module.c"
SETUP_PATH = pathlib.Path(__file__).resolve().parent.parent / "setup.py"


# ------------------------------------------------------------------
# Patch detection and application.
#
# The vendored files carry no local patches as of q3c v2.0.6 (see PATCHES in the
# script), so the machinery is exercised on a synthetic patch with two sites: the
# next reachable upstream bug will need it, and it must still refuse to guess.
# ------------------------------------------------------------------

FIXTURE_PATCH = sync_q3c.Patch(
	name="fixture buf[8]",
	filename="fixture.c",
	upstream_text="    int buf[4];\n",
	patched_text="    int buf[8];\n",
	occurrences=2,
	issue="https://example.invalid/issues/1",
	guard_test="test_fixture_guard")

FIXTURE_UPSTREAM = ("void f(void)\n{\n    int buf[4];\n}\n"
					"void g(void)\n{\n    int buf[4];\n}\n")
FIXTURE_PATCHED = FIXTURE_UPSTREAM.replace("buf[4]", "buf[8]")

def test_no_local_patches_are_carried():
	'''
	Test that the script registers exactly the local divergences CLAUDE.md documents,
	which since the v2.0.6 sync is none. A local patch added to the vendored files
	without being registered here would make the next sync's divergence check stop,
	and one registered without being documented would be invisible to a reader.
	'''
	assert sync_q3c.PATCHES == ()

def test_patch_round_trip():
	'''
	Test that applying a patch to the upstream text patches every site, that reversing
	it reproduces the upstream text byte for byte, and that each side reports the
	right state.
	'''
	assert sync_q3c.patch_state(FIXTURE_PATCH, FIXTURE_UPSTREAM) is sync_q3c.PatchState.NEEDED
	patched, state = sync_q3c.apply_patch(FIXTURE_PATCH, FIXTURE_UPSTREAM)
	assert state is sync_q3c.PatchState.NEEDED
	assert patched == FIXTURE_PATCHED
	assert sync_q3c.patch_state(FIXTURE_PATCH, patched) is sync_q3c.PatchState.ALREADY_UPSTREAM
	assert sync_q3c.unapply_patch(FIXTURE_PATCH, patched) == FIXTURE_UPSTREAM

def test_patch_already_upstream_is_a_no_op():
	'''
	Test that a source already containing the fix is left alone rather than patched twice.
	'''
	unchanged, state = sync_q3c.apply_patch(FIXTURE_PATCH, FIXTURE_PATCHED)
	assert state is sync_q3c.PatchState.ALREADY_UPSTREAM
	assert unchanged == FIXTURE_PATCHED

def test_patch_refuses_an_unrecognized_source():
	'''
	Test that a file in which the patch anchor cannot be found raises rather than being
	silently skipped. This is the case where upstream reworked the surrounding code, or
	fixed the bug with different text, as it did for both v2.0.5 patches in v2.0.6.
	'''
	with pytest.raises(sync_q3c.SyncError) as excinfo:
		sync_q3c.patch_state(FIXTURE_PATCH, "int main(void) { return 0; }\n")
	assert FIXTURE_PATCH.issue in str(excinfo.value)
	assert FIXTURE_PATCH.guard_test in str(excinfo.value)

def test_patch_refuses_a_changed_number_of_sites():
	'''
	Test that a source carrying the wrong number of patch sites raises: upstream adding
	a third copy of a patched loop must not have two of its three sites patched.
	'''
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.patch_state(FIXTURE_PATCH, FIXTURE_UPSTREAM + FIXTURE_PATCH.upstream_text)

def test_patch_refuses_a_partially_patched_source():
	'''
	Test that a source containing both forms - as it would if upstream fixed some but
	not all of the sites - raises rather than being patched the rest of the way.
	'''
	partial = FIXTURE_UPSTREAM.replace("buf[4]", "buf[8]", 1)
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.patch_state(FIXTURE_PATCH, partial)

def test_unapply_refuses_an_unpatched_source():
	'''
	Test that reversing a patch that is not there raises; this is the check that keeps
	the "no undocumented divergence" step honest.
	'''
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.unapply_patch(FIXTURE_PATCH, FIXTURE_UPSTREAM)

def test_patch_rejects_overlapping_texts():
	'''
	Test that a patch whose upstream and patched texts contain one another cannot be
	defined; the occurrence counting that decides the patch state depends on the two
	being distinct.
	'''
	with pytest.raises(ValueError):
		sync_q3c.Patch(name="overlap", filename="fixture.c", upstream_text="int buf;\n",
					   patched_text="int buf;\nint guard;\n", occurrences=1,
					   issue="https://example.invalid/issues/2", guard_test="test_overlap")


# ------------------------------------------------------------------
# Comment stripping
# ------------------------------------------------------------------

def test_strip_c_comments_keeps_string_literals():
	'''
	Test that a "//" inside a string literal is not treated as a comment. The wrapper
	has one (a URL in an error message); a naive stripper would drop the rest of that
	line, and with it any identifier on it.
	'''
	source = 'const char *u = "https://example.org/x"; int keep_me;\n'
	stripped = sync_q3c.strip_c_comments(source)
	assert "keep_me" in stripped
	assert "https://example.org/x" in stripped

def test_strip_c_comments_removes_both_comment_forms():
	'''
	Test that line and block comments are removed and do not glue their neighbors together.
	'''
	stripped = sync_q3c.strip_c_comments("int a; // gone\nint/*gone*/b;\n")
	assert "gone" not in stripped
	assert "intb" not in stripped
	assert "int a;" in stripped

def test_strip_c_comments_rejects_an_unterminated_block_comment():
	'''
	Test that a truncated source raises instead of quietly discarding its tail.
	'''
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.strip_c_comments("int a; /* never closed\n")


# ------------------------------------------------------------------
# The C interface diff
# ------------------------------------------------------------------

def test_interface_index_finds_the_real_declarations():
	'''
	Test that the header index picks up the function prototypes, structs and macros
	the wrapper depends on out of the real common.h.
	'''
	index = sync_q3c.interface_index((VENDOR_DIR / "common.h").read_text())
	for name in ("q3c_ang2ipix", "q3c_radial_query", "q3c_check_sphere_point_in_poly",
				 "init_q3c1", "q3c_prm", "q3c_prm_write", "q3c_ipix_t", "q3c_coord_t",
				 "Q3C_NFULLS", "Q3C_NPARTIALS", "Q3C_MAX_N_POLY_VERTEX", "UNWRAP_RA"):
		assert name in index, f"'{name}' was not indexed from common.h"
	assert "hprm" in index["q3c_radial_query"][0]
	assert "const q3c_ipix_t nside" in index["q3c_prm"][0]

def test_interface_index_records_both_arms_of_an_ifdef():
	'''
	Test that a macro defined twice (common.h defines its math macros once per branch
	of the Q3C_LONG_DOUBLE #ifdef) keeps both definitions, so a change to either is seen.
	'''
	index = sync_q3c.interface_index((VENDOR_DIR / "common.h").read_text())
	assert len(index["q3c_sqrt"]) == 2
	assert index["q3c_sqrt"] != index["q3c_atan"]

def test_wrapper_identifiers_ignore_comments_and_strings():
	'''
	Test that a name mentioned only in a comment or a message string is not counted as
	a use of the C interface.
	'''
	names = sync_q3c.wrapper_identifiers('// q3c_commented\nint x = 1; const char *m = "q3c_quoted";\n')
	assert "x" in names
	assert "q3c_commented" not in names
	assert "q3c_quoted" not in names

def test_no_interface_differences_against_itself():
	'''
	Test that the real header compared with itself reports nothing.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	assert sync_q3c.interface_differences(header, header, WRAPPER_PATH.read_text()) == []

def test_reformatting_alone_is_not_an_interface_difference():
	'''
	Test that rewrapping a prototype across lines is not reported; upstream reformats
	its sources from time to time and a whitespace-only change is not an interface change.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	old = "void q3c_ipix2ang(struct q3c_prm *, q3c_ipix_t, q3c_coord_t *, q3c_coord_t *);"
	assert old in header
	reformatted = header.replace(old, "void q3c_ipix2ang(struct q3c_prm *,\n                  q3c_ipix_t,\n"
									  "                  q3c_coord_t *,\n                  q3c_coord_t *);")
	assert sync_q3c.interface_differences(header, reformatted, WRAPPER_PATH.read_text()) == []

def test_a_changed_prototype_is_reported():
	'''
	Test that adding a parameter to a function the wrapper calls is caught.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	changed = header.replace("void q3c_ipix2ang(struct q3c_prm *, q3c_ipix_t, q3c_coord_t *, q3c_coord_t *);",
							 "void q3c_ipix2ang(struct q3c_prm *, q3c_ipix_t, q3c_coord_t *, q3c_coord_t *, int);")
	assert changed != header
	differences = sync_q3c.interface_differences(header, changed, WRAPPER_PATH.read_text())
	assert len(differences) == 1
	assert "q3c_ipix2ang" in differences[0]

def test_a_changed_buffer_size_macro_is_reported():
	'''
	Test that a change to one of the constants sizing the wrapper's stack buffers is
	caught. Q3C_NFULLS sizes an array in qlsc_q3c_radial_query; a silent change there
	would be a stack overflow.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	changed = header.replace("#define Q3C_NFULLS 50", "#define Q3C_NFULLS 80")
	assert changed != header
	differences = sync_q3c.interface_differences(header, changed, WRAPPER_PATH.read_text())
	assert len(differences) == 1
	assert "Q3C_NFULLS" in differences[0]

def test_a_changed_struct_is_reported():
	'''
	Test that a change to the q3c_prm struct, whose members the wrapper reads and frees
	directly, is caught.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	changed = header.replace("    const q3c_ipix_t *ybits1;\n};\nstruct q3c_prm_write",
							 "    const q3c_ipix_t *ybits1;\n    const q3c_ipix_t *zbits1;\n};\nstruct q3c_prm_write")
	assert changed != header
	differences = sync_q3c.interface_differences(header, changed, WRAPPER_PATH.read_text())
	assert len(differences) == 1
	assert "q3c_prm" in differences[0]

def test_a_removed_declaration_is_reported():
	'''
	Test that a function the wrapper calls disappearing from the header is caught.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	changed = header.replace("q3c_coord_t q3c_sindist(q3c_coord_t, q3c_coord_t, q3c_coord_t, q3c_coord_t);", "")
	assert changed != header
	differences = sync_q3c.interface_differences(header, changed, WRAPPER_PATH.read_text())
	assert len(differences) == 1
	assert "q3c_sindist" in differences[0]
	assert "no longer defined" in differences[0]

def test_changes_the_wrapper_does_not_use_are_ignored():
	'''
	Test that upstream is free to change anything qlsc never calls.
	'''
	header = (VENDOR_DIR / "common.h").read_text()
	assert "q3c_get_nearby_split" not in sync_q3c.wrapper_identifiers(WRAPPER_PATH.read_text())
	changed = header.replace("void q3c_get_nearby_split(struct q3c_prm *, q3c_coord_t, q3c_coord_t, "
							 "q3c_coord_t, q3c_ipix_t *, int);",
							 "void q3c_get_nearby_split(struct q3c_prm *, q3c_coord_t, int);")
	assert changed != header
	assert sync_q3c.interface_differences(header, changed, WRAPPER_PATH.read_text()) == []


# ------------------------------------------------------------------
# Version handling
# ------------------------------------------------------------------

def test_version_from_tag():
	'''
	Test that the Q3C_VERSION string is derived from a tag name.
	'''
	assert sync_q3c.version_from_tag("v2.0.5") == "2.0.5"
	assert sync_q3c.version_from_tag("2.0.6") == "2.0.6"

def test_version_from_an_unparseable_tag_raises():
	'''
	Test that a tag the version cannot be read out of stops the sync rather than
	writing a guess into setup.py.
	'''
	for tag in ("release-2.0.5", "v2.0.5-rc1", "main", ""):
		with pytest.raises(sync_q3c.SyncError):
			sync_q3c.version_from_tag(tag)

def test_setup_q3c_version_matches_the_vendored_sources():
	'''
	Test that the version string can be read out of the real setup.py.
	'''
	assert sync_q3c.setup_q3c_version(SETUP_PATH.read_text()) == "2.0.6"

def test_set_setup_q3c_version_changes_only_that_string():
	'''
	Test that writing a new version rewrites exactly one line of setup.py.
	'''
	original = SETUP_PATH.read_text()
	updated = sync_q3c.set_setup_q3c_version(original, "2.0.7")
	assert sync_q3c.setup_q3c_version(updated) == "2.0.7"
	differing = [(a, b) for a, b in zip(original.splitlines(), updated.splitlines()) if a != b]
	assert len(differing) == 1
	assert "Q3C_VERSION" in differing[0][0]
	assert len(original.splitlines()) == len(updated.splitlines())

def test_setup_version_missing_raises():
	'''
	Test that a setup.py whose define_macros entry has been reworded stops the sync.
	'''
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.setup_q3c_version("define_macros = []\n")

def test_control_version():
	'''
	Test that upstream's own version string is read from a q3c.control file.
	'''
	control = "comment = 'q3c sky indexing plugin'\ndefault_version = '2.0.5'\nrelocatable = true\n"
	assert sync_q3c.control_version(control) == "2.0.5"
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.control_version("comment = 'q3c'\n")
