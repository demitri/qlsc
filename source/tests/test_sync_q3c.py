
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
# Patch detection against the real vendored sources.
#
# Reversing each patch produces the upstream text it was made against, so the real
# files serve as their own fixtures: the "upstream" side is derived from them rather
# than being a hand-copied snippet that could drift out of step with the C code.
# ------------------------------------------------------------------

@pytest.mark.parametrize("patch", sync_q3c.PATCHES, ids=lambda p: p.name)
def test_patch_is_present_in_the_vendored_source(patch):
	'''
	Test that each patch this script knows about is actually in the vendored file,
	i.e. that the script's idea of the local divergences matches the C sources.
	'''
	source = (VENDOR_DIR / patch.filename).read_text()
	assert sync_q3c.patch_state(patch, source) is sync_q3c.PatchState.ALREADY_UPSTREAM
	assert source.count(patch.patched_text) == patch.occurrences

@pytest.mark.parametrize("patch", sync_q3c.PATCHES, ids=lambda p: p.name)
def test_patch_round_trip(patch):
	'''
	Test that reversing a patch and re-applying it reproduces the vendored file byte
	for byte, and that the reversed text is correctly reported as still needing it.
	'''
	vendored = (VENDOR_DIR / patch.filename).read_text()

	upstream_form = sync_q3c.unapply_patch(patch, vendored)
	assert upstream_form != vendored
	assert sync_q3c.patch_state(patch, upstream_form) is sync_q3c.PatchState.NEEDED

	repatched, state = sync_q3c.apply_patch(patch, upstream_form)
	assert state is sync_q3c.PatchState.NEEDED
	assert repatched == vendored

@pytest.mark.parametrize("patch", sync_q3c.PATCHES, ids=lambda p: p.name)
def test_patch_already_upstream_is_a_no_op(patch):
	'''
	Test that a source already containing the fix is left alone rather than patched twice.
	'''
	vendored = (VENDOR_DIR / patch.filename).read_text()
	unchanged, state = sync_q3c.apply_patch(patch, vendored)
	assert state is sync_q3c.PatchState.ALREADY_UPSTREAM
	assert unchanged == vendored

@pytest.mark.parametrize("patch", sync_q3c.PATCHES, ids=lambda p: p.name)
def test_patch_refuses_an_unrecognized_source(patch):
	'''
	Test that a file in which the patch anchor cannot be found raises rather than being
	silently skipped. This is the case where upstream reworked the surrounding code.
	'''
	with pytest.raises(sync_q3c.SyncError) as excinfo:
		sync_q3c.patch_state(patch, "int main(void) { return 0; }\n")
	assert patch.issue in str(excinfo.value)
	assert patch.guard_test in str(excinfo.value)

@pytest.mark.parametrize("patch", sync_q3c.PATCHES, ids=lambda p: p.name)
def test_patch_refuses_a_changed_number_of_sites(patch):
	'''
	Test that a source carrying the wrong number of patch sites raises. Upstream adding
	a fourth copy of the radial query loop must not have three of its four sites patched.
	'''
	upstream_form = sync_q3c.unapply_patch(patch, (VENDOR_DIR / patch.filename).read_text())
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.patch_state(patch, upstream_form + patch.upstream_text)

@pytest.mark.parametrize("patch", sync_q3c.PATCHES, ids=lambda p: p.name)
def test_patch_refuses_a_partially_patched_source(patch):
	'''
	Test that a source containing both forms - as it would if upstream fixed some but
	not all of the sites - raises rather than being patched the rest of the way.
	'''
	vendored = (VENDOR_DIR / patch.filename).read_text()
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.patch_state(patch, vendored + patch.upstream_text)

def test_unapply_refuses_an_unpatched_source():
	'''
	Test that reversing a patch that is not there raises; this is the check that keeps
	the "no undocumented divergence" step honest.
	'''
	patch = sync_q3c.POLY_POINTS_PATCH
	with pytest.raises(sync_q3c.SyncError):
		sync_q3c.unapply_patch(patch, patch.upstream_text)

def test_patches_cover_both_documented_divergences():
	'''
	Test that the script knows about exactly the two divergences CLAUDE.md documents.
	A third local patch added without being registered here would be silently
	overwritten by the next sync.
	'''
	assert {p.filename for p in sync_q3c.PATCHES} == {"q3c_poly.c", "q3cube.c"}
	assert len(sync_q3c.PATCHES) == 2

def test_patched_and_upstream_texts_are_disjoint():
	'''
	Test that no patch's upstream text is a substring of its own patched text; the
	occurrence counting that decides the patch state depends on the two being distinct.
	'''
	for patch in sync_q3c.PATCHES:
		assert patch.upstream_text not in patch.patched_text
		assert patch.patched_text not in patch.upstream_text


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
	assert sync_q3c.setup_q3c_version(SETUP_PATH.read_text()) == "2.0.5"

def test_set_setup_q3c_version_changes_only_that_string():
	'''
	Test that writing a new version rewrites exactly one line of setup.py.
	'''
	original = SETUP_PATH.read_text()
	updated = sync_q3c.set_setup_q3c_version(original, "2.0.6")
	assert sync_q3c.setup_q3c_version(updated) == "2.0.6"
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
