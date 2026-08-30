#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/ndarrayobject.h> // needed to work with NumPy arrays

#include "q3c/common.h"
//#include "my_bits.h"

// character parameter values
// ref: https://docs.python.org/3/c-api/arg.html?highlight=pyarg_parsetupleandkeywords#building-values

#define Q3C_STRUCT_POINTER_BUFFER "Q3C_prm_struct_pointer"

// The largest semi-major axis an elliptical region may have; see the comment in
// validate_ellipse_parameters for why. Exported to Python as
// qlsc.q3c.ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS so the two layers cannot drift apart.
#define ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS 90.0

// ------------------------------------------------------------------
// Right ascension conventions in the vendored Q3C code
// ------------------------------------------------------------------
//
// Several vendored functions carry the documented precondition "strictly
// 0<=ra<360" and do nothing to enforce it: q3c_get_facenum computes the cube
// face as (char)fmod((ra + 45) / 90, 4) + 1, which for a negative ra truncates
// toward zero and can return a face number that is not even in [0,5] (ra=-315
// returns -2). Everything downstream of a wrong face is wrong or spuriously
// rejected as "too large".
//
// This wrapper is the boundary where that is fixed - the vendored files are not
// edited for it. Every entry point below that hands an ra to Q3C is audited here;
// re-run this audit when a Q3C function is newly exported or the vendored code is
// synced (scripts/sync_q3c.py diffs the declarations this depends on).
//
//   SAFE, unwraps internally:
//     q3c_ang2ipix, q3c_ang2ipix_xy  - both call UNWRAP_RA(ra0) on entry
//   SAFE, periodic in ra with period 360:
//     q3c_dist, q3c_sindist  - ra enters only as sin((ra1-ra2)/2)^2
//     q3c_in_ellipse         - ra enters only as cos/sin of (alpha1-alpha)
//     q3c_xy2ang, q3c_xy2facenum, q3c_ipix2ang, q3c_ipix2xy, q3c_pixarea
//                            - take no ra at all
//   NEEDS UNWRAPPING BY THIS WRAPPER:
//     q3c_get_facenum        - unwrapped in qlsc_q3c_facenum
//     q3c_radial_query       - unwrapped in qlsc_q3c_radial_query
//     q3c_ellipse_query      - unwrapped in qlsc_q3c_ellipse_query
//     q3c_check_sphere_point_in_poly, q3c_poly_query
//                            - both reach q3c_get_facenum through
//                              q3c_get_facenum_poly, which reads the FIRST
//                              polygon vertex; the vertex ra values are unwrapped
//                              into a local buffer by polygon_vertex_views

void *get_item_pointer(int ndim, void *buf, Py_ssize_t *strides,
   Py_ssize_t *suboffsets, Py_ssize_t *indices) {
   char *pointer = (char*)buf;
   int i;
   for (i = 0; i < ndim; i++) {
	   pointer += strides[i] * indices[i];
	   if (suboffsets[i] >=0 ) {
		   pointer = *((char**)pointer) + suboffsets[i];
	   }
   }
   return (void*)pointer;
}

// ---------------------------------
// prm (main Q3C structure) routines
// ---------------------------------
// TODO: submit this as a pull request to include in q3cube.c ?
void free_q3c1(struct q3c_prm *hprm)
{
	// the struct members are declared const in q3c v2.0.1+; they were malloc'd via init_q3c1
	free((void *)hprm->xbits);
	free((void *)hprm->ybits);
	free((void *)hprm->xbits1);
	free((void *)hprm->ybits1);
	free(hprm);
}

// Function that Python will call when the prm (wrapped in a PyObject) is freed.
static void del_prm(PyObject *obj)
{
	// get C pointer out of object
	struct q3c_prm *hprm = PyCapsule_GetPointer(obj, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL) {
		// cannot raise from a destructor; report the unraisable error instead
		PyErr_WriteUnraisable(obj);
		return;
	}
	free_q3c1(hprm);
}

static PyObject *
qlsc_init_q3c1(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *returnValue = NULL;
	
	// dynamically allocate memory for this as it will be kept around
	// (q3c v2.0.1+ initializes into a q3c_prm_write, which has the same
	// layout as q3c_prm without the const qualifiers)
	struct q3c_prm_write *hprm = malloc(sizeof (struct q3c_prm_write));
	q3c_ipix_t nside = 0; // no default value
	
	// keyword list
	static char *kwlist[] = {"nside", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "L", 			// "L" = long long
									 kwlist,
									 &nside))
	{
		goto except;
	}
	init_q3c1(hprm, nside);
	
	//PySys_WriteStdout("hprm = %"PRId64"\n", hprm->nside);
	
	// return C pointer wrapped in Python object
	returnValue = PyCapsule_New(hprm, Q3C_STRUCT_POINTER_BUFFER, del_prm); // provide function used to free pointer
	goto finally;
	
except:
	returnValue = NULL;
finally:
	return returnValue;
}

// ---------------------------------

static PyObject *
qlsc_q3c_nside(PyObject *self, PyObject *arg)
{
	PyObject *return_value = NULL;

	struct q3c_prm *hprm;
	
	// no parsing needed; one value is expected
	//hprm = (struct q3c_prm*)PyCapsule_GetPointer(arg, Q3C_STRUCT_POINTER_BUFFER);
	
	PyObject *capsule;
	if (!PyArg_ParseTuple(arg, "O", &capsule))
		return NULL;
	hprm = PyCapsule_GetPointer(capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	return_value = PyLong_FromLongLong(hprm->nside);
	
//	goto finally;
//except:
//	return_value = NULL;
//finally:
	return return_value;
}

static PyObject *
qlsc_q3c_ang2ipix(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
//qlsc_q3c_ang2ipix(PyObject *module, PyObject *args) // -> cast as PyCFunction
{
	// external parameters
	// -------------------
	PyObject *hprm_capsule;
	q3c_coord_t ra, dec;
	
	// used for numpy values
	// ---------------------
	PyObject *np_ra, *np_dec;

	// internal variables
	// ------------------
	struct q3c_prm *hprm;
	int array_form = 0; // set to 1 if the inputs are numpy arrays
	q3c_ipix_t ipix = 0;
	//Py_ssize_t n; // number of points

	// returned values
	// ---------------
	PyArrayObject *np_ipix;
	
	static char *kwlist[] = {"hprm", "ra", "dec", NULL};

	// We accept one of two parameter signatures:
	//  ang2ipix ( <object>, <double>, <double> )
	// and
	//  ang2ipix ( <object>, <ndarray:double>, <ndarray:double> )

	//PySys_WriteStdout("about to parse input\n");
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "Odd",	// object + 2 doubles
    								 kwlist,
    								 &hprm_capsule, &ra, &dec))
	{
		// AN ERROR HAS BEEN SET HERE!
		// Clear it since we're going to try to continue.
		PyErr_Clear();
		
		// unable to parse inputs
		// not two doubles; see if we've been given NumPy arrays:
		if (PyArg_ParseTupleAndKeywords(args, kwargs,
										 "OOO",	// 3 objects
										 kwlist,
										 &hprm_capsule, &np_ra, &np_dec))
		{
			// success - found three objects
			array_form = 1;
		} else {
			PyErr_SetString(PyExc_TypeError, "Could not determine type for ra or dec: either use scalar values or NumPy arrays.");
			return NULL; // -> raise exception
		}
	}

	// TODO: look for examples where attribute names are pre-cached.
	//PyObject *hprm_attr_name = PyUnicode_FromString("_hprm");
	//PySys_WriteStdout("has attr: %d\n", PyObject_HasAttr(self, hprm_attr_name));

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	if (array_form) {
		
		Py_buffer ra_view, dec_view;
		q3c_ipix_t *ipix_array; // hold the calculated ipix values

		// NUMPY ARRAY VERSION
		
		// Good reference for this stuff: https://documentation.help/Python-3.7/buffer.html
		
		//
		// Create views from the numpy objects from which we can access the C array and other info.
		// Ensure that the memory is contigguous and that the object reports its data type.
		// Successful calls to PyObject_GetBuffer must be balanced by PyBuffer_Release!
		// PyObject_GetBuffer returns 0 on success, -1 otherwise.
		//
		// PyBUF_FULL_RO = (PyBUF_INDIRECT | PyBUF_STRIDES | PyBUF_FORMAT)
		//
		if (PyObject_GetBuffer(np_ra, &ra_view, PyBUF_FULL_RO) == -1)
			return NULL;
		if (PyObject_GetBuffer(np_dec, &dec_view, PyBUF_FULL_RO) == -1) {
			PyBuffer_Release(&ra_view);
			return NULL;
		}

		// check the data type is double for ra,dec
		if (strcmp(ra_view.format, "d") != 0 || strcmp(dec_view.format, "d") != 0) {
			PyErr_SetString(PyExc_TypeError, "'ra' and 'dec' arrays must be of type np.double.");
			goto ang2ipix_array_early_exit;
		}
		
		// check arrays are the same length
		// Py_ssize_t *shape
		if (ra_view.shape[0] != dec_view.shape[0]) {
			PyErr_SetString(PyExc_ValueError, "'ra' and 'dec' arrays must be the same length.");
			goto ang2ipix_array_early_exit;
		}

		// check that the arrays are 1D - don't do this here; we'll accept slices where the array is complex
		//		if ((ra_view.ndim != 1) || (ra_view.ndim != 1)) {
		//			PyErr_SetString(PyExc_ValueError, "'ra' and 'dec' must be 1D arays.");
		//			goto ang2ipix_array_early_exit;
		//		}

		Py_ssize_t n = ra_view.shape[0]; // number of points

		// Debugging.
//		Py_ssize_t *strides = ra_view.strides;
//		PySys_WriteStdout("ndim: %d\n", ra_view.ndim);
//		PySys_WriteStdout("ra,dec: dim[0] = (%zd,%zd)\n", ra_view.shape[0], dec_view.shape[0]);
//		PySys_WriteStdout("stride: %zd\n", strides[0]);
//		PySys_WriteStdout("contiguous? %d\n", strides[0] == sizeof(double));
//		if (ra_view.suboffsets)
//			PySys_WriteStdout("suboffset: %zd\n", ra_view.suboffsets[0]);
		
//		double *arr = ra_view.buf;
//		PySys_WriteStdout("First element %f\n", arr[0]);
//		for (int i=0; i < n; i++) {
//			PySys_WriteStdout("%f\n", *(double*)((char*)arr + i * strides[0]));
//		}
		
		q3c_coord_t *ra_buffer = ra_view.buf; // pointer to C array
		q3c_coord_t *dec_buffer = dec_view.buf;

		if (ra_view.suboffsets || dec_view.suboffsets) {
			// complex array! requires special handling
			
			PyErr_SetString(PyExc_NotImplementedError, "This is a particularly complex array! Please send reproducible code to the author to support this type of array. In the meantime, pass a simpler (e.g. 1D) array to this function.");
			goto ang2ipix_array_early_exit;

			// The documentation says these are PIL-style arrays (Python Imaging Library).
			// This seems way out of scope for this kind of function, but the code
			// would look something like this. Not clear what the "indices" paramter is.
			//
			// An example of usage is here: https://github.com/python/cpython/blob/master/Objects/abstract.c#L611
			//
			// This is a start on the code, but it's more complicated than this.
//			q3c_coord_t *ra_ptr;
//			q3c_coord_t *dec_ptr;
//			for (int i=0; i < n; i++) {
//				ra_ptr  = get_item_pointer(ra_view.ndim, ra_buffer, ra_view.strides, ra_view.suboffsets, <indices>??);
//				dec_ptr = get_item_pointer(dec_view.ndim, dec_buffer, dec_view.strides, dec_view.suboffsets, <indices>??);
//				q3c_ang2ipix(hprm, *((double *)ra_ptr), *((double *)dec_ptr), &ipix_array[i]);
//			}
			
		} else if ((ra_view.strides[0] != sizeof(q3c_coord_t)) || dec_view.strides[0] != sizeof(q3c_coord_t)) {
			// handle strides
			
			ipix_array = malloc(n * sizeof(q3c_ipix_t)); // allocate the ipix values array
			
			Py_ssize_t ra_stride  = ra_view.strides[0];
			Py_ssize_t dec_stride = dec_view.strides[0];
			q3c_coord_t ra, dec;
			for (int i=0; i < n; i++) {
				ra = *(double*)((char*)ra_buffer + i * ra_stride);
				dec = *(double*)((char*)dec_buffer + i * dec_stride);
				// non-finite values reach an undefined float->int conversion (upstream q3c issue #52)
				if (!isfinite(ra) || !isfinite(dec)) {
					PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' arrays must contain only finite values.");
					free(ipix_array);
					goto ang2ipix_array_early_exit;
				}
				q3c_ang2ipix(hprm, ra, dec, &ipix_array[i]);
			}

		} else {
			// standard C array, contiguous data

			ipix_array = malloc(n * sizeof(q3c_ipix_t)); // allocate the ipix values array

			for (int i=0; i < n; i++) {
				// non-finite values reach an undefined float->int conversion (upstream q3c issue #52)
				if (!isfinite(ra_buffer[i]) || !isfinite(dec_buffer[i])) {
					PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' arrays must contain only finite values.");
					free(ipix_array);
					goto ang2ipix_array_early_exit;
				}
				q3c_ang2ipix(hprm, ra_buffer[i], dec_buffer[i], &ipix_array[i]);
			}
		}

		// create a new NumPy array with the ipix values
		npy_intp dims[1];
		dims[0] = n;
		np_ipix = (PyArrayObject *)PyArray_SimpleNewFromData(1, 			// int nd - length of dims array
															 dims, 			// npy_intp* dims (for 1D, could just be an int)
															 NPY_INT64, 	// int typenum,
															 ipix_array); 	// void* data
		PyArray_ENABLEFLAGS(np_ipix, NPY_ARRAY_OWNDATA);

		PyBuffer_Release(&ra_view);
		PyBuffer_Release(&dec_view);
		
		// if there is an error at this point, print it
		//PyErr_Print();
		
		// return numpy array ("N" steals the reference; "O" would leak one per call)
		return Py_BuildValue("N", np_ipix);

ang2ipix_array_early_exit:
		// release references, bail
		PyBuffer_Release(&ra_view);
		PyBuffer_Release(&dec_view);
		return NULL;
	}
	else {
		// SCALAR VERSION

		// Non-finite values would hit an undefined float->int conversion inside
		// q3c_ang2ipix (upstream q3c issue #52); returning None here (the old
		// behavior) silently hid the bad input.
		if ((!isfinite(ra)) || (!isfinite(dec)))
		{
			PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' values must be finite.");
			return NULL;
		}

		q3c_ang2ipix(hprm, ra, dec, &ipix);

		return PyLong_FromLongLong(ipix);
	}

}

static PyObject *
qlsc_q3c_ang2ipix_xy(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
{
	// external parameters
	q3c_coord_t ra;
	q3c_coord_t dec;
	PyObject *hprm_capsule;

	//internal variables
	struct q3c_prm *hprm;

	// return values
	char facenum;
	q3c_ipix_t ipix = 0;
	q3c_coord_t x, y;
	
	static char *kwlist[] = {"hprm", "ra", "dec", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "Odd",	// object + 2 doubles
    								 kwlist,
    								 &hprm_capsule, &ra, &dec))
	{
		// unable to parse inputs -> raise exception
        return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	// Non-finite values would hit an undefined float->int conversion inside
	// q3c_ang2ipix_xy (upstream q3c issue #52); returning None here (the old
	// behavior) silently hid the bad input.
	if ((!isfinite(ra)) || (!isfinite(dec)))
	{
		PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' values must be finite.");
		return NULL;
	}

	// q3c_ang2ipix_xy (struct q3c_prm *hprm, q3c_coord_t ra0, q3c_coord_t dec0,
	//					char *out_face_num, q3c_ipix_t *ipix, q3c_coord_t *x_out,
	//					q3c_coord_t *y_out)
	q3c_ang2ipix_xy (hprm, ra, dec,
					 &facenum,
					 &ipix,
					 &x, &y);
	
	// return values as dictionary
	PyObject* return_dict = PyDict_New();

	PyObject* facenum_py = PyLong_FromLong((long)facenum);
	PyObject* ipix_py = PyLong_FromLongLong(ipix);
	PyObject* x_py = PyFloat_FromDouble(x);
	PyObject* y_py = PyFloat_FromDouble(y);

	PyDict_SetItemString(return_dict, "facenum", facenum_py);
	PyDict_SetItemString(return_dict, "x", x_py);
	PyDict_SetItemString(return_dict, "y", y_py);
	PyDict_SetItemString(return_dict, "ipix", ipix_py);

	Py_DECREF(facenum_py);
	Py_DECREF(ipix_py);
	Py_DECREF(x_py);
	Py_DECREF(y_py);

	return return_dict;
}

// ------------------------------------------------------------------
// Shared helpers for the ipix2* functions, which accept either a
// single ipix value or an array of them.
// ------------------------------------------------------------------

// Validate a single ipix value against the scheme described by hprm.
// Returns 0 if the value is in range, -1 with an exception set otherwise.
// Out-of-range values are not silently clamped: q3c_ipix2ang would happily
// decode them into coordinates that belong to no pixel of this scheme.
static int
ipix_bounds_check(struct q3c_prm *hprm, q3c_ipix_t ipix)
{
	// nside <= 2^30, so 6 * nside^2 <= 6.92e18 and cannot overflow int64
	q3c_ipix_t nbins = 6 * hprm->nside * hprm->nside;
	if (ipix < 0 || ipix >= nbins) {
		PyErr_Format(PyExc_ValueError,
					 "The ipix value %lld is out of range for this scheme; expected [0,%lld].",
					 (long long)ipix, (long long)(nbins - 1));
		return -1;
	}
	return 0;
}

// True if the given object should be handled as an array of ipix values
// rather than as a single value. A 0-d NumPy array is a scalar here, matching
// the convention ang2ipix uses for its ra,dec arguments.
static int
is_ipix_array_input(PyObject *obj)
{
	if (PyArray_Check(obj))
		return PyArray_NDIM((PyArrayObject *)obj) > 0;
	if (PyList_Check(obj) || PyTuple_Check(obj))
		return 1;
	return 0;
}

// Extract a single ipix value from a Python object.
// Returns 0 on success, -1 with an exception set otherwise. PyNumber_Index
// accepts Python and NumPy integers (and 0-d integer arrays) but rejects
// floats, so a fractional pixel number fails loudly rather than truncating.
static int
ipix_from_object(PyObject *obj, q3c_ipix_t *ipix)
{
	PyObject *index = PyNumber_Index(obj);
	if (index == NULL)
		return -1;
	long long value = PyLong_AsLongLong(index);
	Py_DECREF(index);
	if (value == -1 && PyErr_Occurred())
		return -1;
	*ipix = (q3c_ipix_t)value;
	return 0;
}

// Convert an array-like object to a C-contiguous, aligned int64 array of ipix
// values, validating every element against the scheme. Returns a new reference
// on success, NULL with an exception set otherwise.
//
// The dtype of the input is settled BEFORE anything is converted. Asking NumPy
// for an int64 array directly is not enough: the "safe cast" rule that refuses
// an existing float64 array does not apply to a value NumPy is constructing for
// the first time, so [1.5, 2.5], [True, False] and ["10", "20"] would all be
// built into an int64 array without a word - 1.5 silently becoming pixel 1.
// Inferring the dtype first and requiring an integer type rejects all three.
static PyArrayObject *
ipix_array_from_object(struct q3c_prm *hprm, PyObject *obj)
{
	// dtype NULL: infer it from the object rather than imposing one
	PyArrayObject *source = (PyArrayObject *)PyArray_FromAny(obj, NULL, 0, 0, 0, NULL);
	if (source == NULL)
		return NULL;

	// PyTypeNum_ISINTEGER spans the signed and unsigned integer types and
	// excludes bool, float, and everything else
	if (!PyTypeNum_ISINTEGER(PyArray_TYPE(source))) {
		PyErr_Format(PyExc_TypeError,
					 "The ipix values must be integers; was given an array of dtype '%S'.",
					 (PyObject *)PyArray_DESCR(source));
		Py_DECREF(source);
		return NULL;
	}

	// NPY_ARRAY_IN_ARRAY handles strides and non-contiguous input; the source is
	// known to be an integer array by this point, so the only conversion that is
	// not value-preserving is uint64 above INT64_MAX, and those values land in
	// the negative half of int64 where the bounds check below rejects them.
	PyArrayObject *array = (PyArrayObject *)PyArray_FROM_OTF((PyObject *)source, NPY_INT64,
															NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST);
	Py_DECREF(source);
	if (array == NULL)
		return NULL;

	q3c_ipix_t *data = (q3c_ipix_t *)PyArray_DATA(array);
	npy_intp n = PyArray_SIZE(array);
	for (npy_intp i = 0; i < n; i++) {
		if (ipix_bounds_check(hprm, data[i]) != 0) {
			Py_DECREF(array);
			return NULL;
		}
	}
	return array;
}

static PyObject *
qlsc_q3c_ipix2ang(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
{
	// external parameters
	PyObject *hprm_capsule;
	PyObject *ipix_obj;

	// internal variables
	struct q3c_prm *hprm;
	q3c_ipix_t ipix; // int64_t
	q3c_coord_t ra;  // double
	q3c_coord_t dec;

	static char *kwlist[] = {"hprm", "ipix", NULL};

	// We accept one of two parameter signatures:
	//  ipix2ang ( <object>, <int> )        -> (ra,dec) as floats
	//  ipix2ang ( <object>, <array-like> ) -> (ra,dec) as two arrays

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "OO",			// object + object
    								 kwlist,
    								 &hprm_capsule, &ipix_obj))
	{
		// unable to parse inputs -> raise exception
        return NULL;
	}

	//PyObject *hprm_attr_name = PyUnicode_FromString("_hprm");
	//PySys_WriteStdout("has attr: %d\n", PyObject_HasAttr(self, hprm_attr_name));

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	if (is_ipix_array_input(ipix_obj)) {
		// ARRAY VERSION

		PyArrayObject *ipix_array = ipix_array_from_object(hprm, ipix_obj);
		if (ipix_array == NULL)
			return NULL;

		// the outputs mirror the shape of the input; NumPy owns this memory
		PyObject *np_ra = PyArray_SimpleNew(PyArray_NDIM(ipix_array), PyArray_DIMS(ipix_array), NPY_DOUBLE);
		if (np_ra == NULL) {
			Py_DECREF(ipix_array);
			return NULL;
		}
		PyObject *np_dec = PyArray_SimpleNew(PyArray_NDIM(ipix_array), PyArray_DIMS(ipix_array), NPY_DOUBLE);
		if (np_dec == NULL) {
			Py_DECREF(ipix_array);
			Py_DECREF(np_ra);
			return NULL;
		}

		q3c_ipix_t *ipix_data = (q3c_ipix_t *)PyArray_DATA(ipix_array);
		q3c_coord_t *ra_data  = (q3c_coord_t *)PyArray_DATA((PyArrayObject *)np_ra);
		q3c_coord_t *dec_data = (q3c_coord_t *)PyArray_DATA((PyArrayObject *)np_dec);
		npy_intp n = PyArray_SIZE(ipix_array);

		for (npy_intp i = 0; i < n; i++)
			q3c_ipix2ang(hprm, ipix_data[i], &ra_data[i], &dec_data[i]);

		Py_DECREF(ipix_array);

		// "N" steals the references; "O" would leak one per call
		return Py_BuildValue("NN", np_ra, np_dec);
	}

	// SCALAR VERSION

	if (ipix_from_object(ipix_obj, &ipix) != 0)
		return NULL;
	if (ipix_bounds_check(hprm, ipix) != 0)
		return NULL;

	q3c_ipix2ang(hprm, ipix, &ra, &dec);

	// debug
	//PySys_WriteStdout("ipix=%lld, ra,dec=(%.4f, %.4f)\n", ipix, ra, dec);

	return Py_BuildValue("dd", ra, dec);
}

void ipix_to_xy(struct q3c_prm *hprm, q3c_ipix_t ipix, q3c_coord_t *x, q3c_coord_t *y)
{
	// taken from q3c_ipix2ang in q3cube.c
	//
	q3c_ipix_t ipix1;
	q3c_ipix_t i2, i3, x0, y0;
	// the bit tables are declared const in the q3c_prm struct (q3c v2.0.1+);
	// they are read-only here, so the pointers keep the qualifier
	const q3c_ipix_t *xbits1;
	const q3c_ipix_t *ybits1;
	const q3c_ipix_t ii1 = 1 << (Q3C_INTERLEAVED_NBITS / 2);

	q3c_ipix_t nside = hprm->nside;
	xbits1 = hprm->xbits1;
	ybits1 = hprm->ybits1;

	ipix1 = ipix % (nside * nside);

	i3 = ipix1 % Q3C_I1;
	i2 = ipix1 / Q3C_I1;
	x0 = xbits1[i3];
	y0 = ybits1[i3];
	i3 = i2 % Q3C_I1;
	i2 = i2 / Q3C_I1;
	x0 += xbits1[i3] * ii1;
	y0 += ybits1[i3] * ii1;
	i3 = i2 % Q3C_I1;
	i2 = i2 / Q3C_I1;
	x0 += xbits1[i3] * ii1 * ii1;
	y0 += ybits1[i3] * ii1 * ii1;
	i3 = i2 % Q3C_I1;
	i2 = i2 / Q3C_I1;
	x0 += xbits1[i3] * ii1 * ii1 * ii1;
	y0 += ybits1[i3] * ii1 * ii1 * ii1;
	
	*x = (((q3c_coord_t)x0) / nside) * 2 - 1;
	*y = (((q3c_coord_t)y0) / nside) * 2 - 1;
	/* Now -1<x<1 and -1<y<1 */
}

static PyObject *
qlsc_q3c_ipix2xy(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	PyObject *ipix_obj;

	// internal variables
	struct q3c_prm *hprm;
	q3c_ipix_t ipix;
//	const q3c_ipix_t ii1 = 1 << (Q3C_INTERLEAVED_NBITS / 2);
//	q3c_ipix_t nside;
//	q3c_ipix_t ipix1;
//	q3c_ipix_t *xbits1;
//	q3c_ipix_t *ybits1;
//	q3c_ipix_t i2, i3, x0, y0;

	// return values
	//char face_num;
	int face_num;
	q3c_coord_t x, y;

	static char *kwlist[] = {"hprm", "ipix", NULL};

	// We accept one of two parameter signatures:
	//  ipix2xy ( <object>, <int> )        -> (facenum,x,y) as scalars
	//  ipix2xy ( <object>, <array-like> ) -> (facenum,x,y) as three arrays

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "OO",						// object + object
    								 kwlist,					//
    								 &hprm_capsule, &ipix_obj))
	{
		// unable to parse inputs -> raise exception
        return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	if (is_ipix_array_input(ipix_obj)) {
		// ARRAY VERSION

		PyArrayObject *ipix_array = ipix_array_from_object(hprm, ipix_obj);
		if (ipix_array == NULL)
			return NULL;

		int ndim = PyArray_NDIM(ipix_array);
		npy_intp *dims = PyArray_DIMS(ipix_array);

		// the outputs mirror the shape of the input; NumPy owns this memory
		PyObject *np_facenum = PyArray_SimpleNew(ndim, dims, NPY_INT64);
		PyObject *np_x = NULL;
		PyObject *np_y = NULL;
		if (np_facenum != NULL)
			np_x = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
		if (np_x != NULL)
			np_y = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
		if (np_y == NULL) {
			Py_DECREF(ipix_array);
			Py_XDECREF(np_facenum);
			Py_XDECREF(np_x);
			return NULL;
		}

		q3c_ipix_t *ipix_data    = (q3c_ipix_t *)PyArray_DATA(ipix_array);
		q3c_ipix_t *facenum_data = (q3c_ipix_t *)PyArray_DATA((PyArrayObject *)np_facenum);
		q3c_coord_t *x_data      = (q3c_coord_t *)PyArray_DATA((PyArrayObject *)np_x);
		q3c_coord_t *y_data      = (q3c_coord_t *)PyArray_DATA((PyArrayObject *)np_y);
		npy_intp n = PyArray_SIZE(ipix_array);

		for (npy_intp i = 0; i < n; i++) {
			// integer arithmetic - see the note in the scalar branch below
			facenum_data[i] = ipix_data[i] / (hprm->nside * hprm->nside);
			ipix_to_xy(hprm, ipix_data[i], &x_data[i], &y_data[i]);
		}

		Py_DECREF(ipix_array);

		// "N" steals the references; "O" would leak one per call
		return Py_BuildValue("NNN", np_facenum, np_x, np_y);
	}

	// SCALAR VERSION

	if (ipix_from_object(ipix_obj, &ipix) != 0)
		return NULL;
	if (ipix_bounds_check(hprm, ipix) != 0)
		return NULL;

	// integer arithmetic - going through doubles (e.g. pow()) rounds ipix values
	// near face boundaries onto the wrong face at depths ≳26
	face_num = (int)(ipix / (hprm->nside * hprm->nside));

//	PySys_WriteStdout("x = %f\n", (double) (ipix / ((hprm->nside)^2)) );
//	PySys_WriteStdout("ipix = %"PRId64"\n", ipix);
//	PySys_WriteStdout("face_num = %d\n", face_num);
//	PySys_WriteStdout("hprm->nside = %f\n", pow(hprm->nside,2));

	ipix_to_xy(hprm, ipix, &x, &y);

	return Py_BuildValue("ldd", face_num, x, y);
}

static PyObject *
qlsc_q3c_pixarea(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external xy2angparameters
	q3c_ipix_t ipix;
	PyObject *hprm_capsule;
	int depth;

	// internal variables
	struct q3c_prm *hprm;
	q3c_coord_t area;
	
	static char *kwlist[] = {"hprm", "ipix", "depth", NULL};
	
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "OLi",						// accept parameters of type double
    								 kwlist,					//
    								 &hprm_capsule, &ipix, &depth))
	{
		// unable to parse inputs -> raise exception
        return NULL;
	}
	
	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception
	
	area = q3c_pixarea(hprm, ipix, depth);

	// debug
	//PySys_WriteStdout("ipix = %"PRId64"\n", ipix);
	//PySys_WriteStdout("depth = %d\n", depth);
	//PySys_WriteStdout("area = %.10e\n", area);
	
	return PyFloat_FromDouble(area);
}

static PyObject *
qlsc_q3c_facenum(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t ra, dec;
	PyObject *hprm_capsule;
	
	// internal variables
	//struct q3c_prm *hprm;
	char facenum;
	
	static char *kwlist[] = {"hprm", "ra", "dec", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									  "Odd",			// object + 2 doubles
									  kwlist,
									  &hprm_capsule, &ra, &dec))
	 {
		 // unable to parse inputs -> raise exception
		 return NULL;
	 }
	
	//hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);

	// Non-finite values reach an undefined float->int conversion in the (char)
	// cast inside q3c_get_facenum (upstream q3c issue #52).
	if ((!isfinite(ra)) || (!isfinite(dec)))
	{
		PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' values must be finite.");
		return NULL;
	}

	// q3c_get_facenum is documented "strictly 0<=ra<360" and does not enforce it;
	// see the audit at the top of this file. Without this, ra=-315 returns -2.
	ra = UNWRAP_RA(ra);

	facenum = q3c_get_facenum(ra, dec);

	return PyLong_FromLong((long)facenum);
}

static PyObject *
qlsc_q3c_dist(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t ra1, dec1; // degrees
	q3c_coord_t ra2, dec2; // degrees
	
	// internal parameters
	q3c_coord_t distance;
	
	static char *kwlist[] = {"ra1", "dec1", "ra2", "dec2", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "dddd",	// accept 4 parameters of type double
									 kwlist,
									 &ra1, &dec1, &ra2, &dec2))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	distance = q3c_dist( ra1, dec1, ra2, dec2 );
	
	return PyFloat_FromDouble(distance);
}

static PyObject *
qlsc_q3c_sindist(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t ra1, dec1; // degrees
	q3c_coord_t ra2, dec2; // degrees
	
	// internal parameters
	q3c_coord_t sine_distance;
	
	static char *kwlist[] = {"ra1", "dec1", "ra2", "dec2", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "dddd",	// accept 4 parameters of type double
									 kwlist,
									 &ra1, &dec1, &ra2, &dec2))
	{
		// unable to parse inputs; PyArg_ParseTupleAndKeywords has set the real
		// exception - do not overwrite it with a generic one
		return NULL;
	}

	sine_distance = q3c_sindist( ra1, dec1, ra2, dec2 );
	
	return PyFloat_FromDouble(sine_distance);
}

static PyObject *
qlsc_q3c_xy2ang(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t x, y;
	char face_num0;

	// returned values
	q3c_coord_t ra = 0;
	q3c_coord_t dec = 0;

	static char *kwlist[] = {"x", "y", "facenum", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "ddb",	// double, double, char
									 kwlist,
									 &x, &y, &face_num0))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	/* This code has been cut out from ipix2ang BEGIN */
	if ((face_num0 >= 1) && (face_num0 <= 4))
	{
		ra = q3c_atan(x);
		dec = Q3C_RADEG * q3c_atan(y * q3c_cos(ra));
		ra = ra * Q3C_RADEG + ((q3c_coord_t)face_num0 - 1) * 90;
		if (ra < 0)
		{
			ra += (q3c_coord_t)360;
		}
	}
	else
	{
		if (face_num0 == 0)
		{
			ra = Q3C_RADEG * q3c_atan2(x, -y);
			dec = Q3C_RADEG * q3c_atan(1 / q3c_sqrt(x * x + y * y));
			if (ra < 0)
			{
				ra += (q3c_coord_t)360;
			}
		}
		if (face_num0 == 5)
		{
			ra = Q3C_RADEG * q3c_atan2(x, y);
			dec = -Q3C_RADEG * q3c_atan(1 / q3c_sqrt(x * x + y * y));
			if (ra < 0)
			{
				ra += (q3c_coord_t)360;
			}

		}
	}
	/* This code has been cut out from ipix2ang END */
	
	return Py_BuildValue("dd", ra, dec);
}

static PyObject *
qlsc_q3c_xy2facenum(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t x, y;
	char facenum;

	static char *kwlist[] = {"x", "y", "facenum", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "ddb",	// double, double, char
									 kwlist,
									 &x, &y, &facenum))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	char result = q3c_xy2facenum(x, y, facenum);

	return PyLong_FromLong( (long) result );
}

// The exception type raised on a q3c internal range-capacity overflow; created in PyInit_q3c.
static PyObject *QlscRangeOverflowError = NULL;

// Raise a RangeOverflowError for an ipix range overflow in one of the q3c region
// query functions. 'query' describes the query as it was actually run (i.e. after
// ra normalization) so that the report is reproducible.
static void raise_query_overflow(const char *function, const char *query)
{
	char msg[1024];
	snprintf(msg, sizeof msg,
			 "%s: too many ipix ranges for the query %s. This is an internal error in the embedded "
			 "Q3C code (version " Q3C_VERSION "); please report it with these values at "
			 "https://github.com/segasai/q3c/issues .",
			 function, query);
	PyErr_SetString(QlscRangeOverflowError, msg);
}

// ------------------------------------------------------------------
// Shared helpers for the region query functions (radial, polygon, ellipse).
// All three hand q3c a pair of fixed-size output arrays and get back ipix
// ranges in the same padded form.
// ------------------------------------------------------------------

// Fill the padded output arrays with q3c's own "no more ranges" filler, so that
// a q3c path returning early without calling array_filler (q3c_poly_query does
// this when the polygon is too large) cannot leave the caller reading
// uninitialized stack memory.
static void
init_padded_ranges(q3c_ipix_t *fulls_padded, q3c_ipix_t *partials_padded)
{
	for (int i = 0; i < 2 * Q3C_NFULLS; i += 2) {
		fulls_padded[i] = 1;
		fulls_padded[i+1] = -1;
	}
	for (int i = 0; i < 2 * Q3C_NPARTIALS; i += 2) {
		partials_padded[i] = 1;
		partials_padded[i+1] = -1;
	}
}

// Build the (fulls, partials) NumPy arrays from the padded arrays q3c fills.
// q3c's array_filler() pads both out to their full length with (1,-1) pairs,
// which are not wanted here; only the exact (1,-1) pair is filler, since a
// legitimate range can start or end at 1.
// Returns a new reference to a 2-tuple of (n,2) int64 arrays, or NULL with an
// exception set.
static PyObject *
ipix_range_arrays(const q3c_ipix_t *fulls_padded, const q3c_ipix_t *partials_padded)
{
	int fulls_filler_pos = 0;
	int partials_filler_pos = 0;
	q3c_ipix_t *fulls, *partials;
	npy_intp dims[2];
	PyArrayObject *np_fulls, *np_partials;

	while (fulls_filler_pos < 2 * Q3C_NFULLS &&
		   !(fulls_padded[fulls_filler_pos] == 1 && fulls_padded[fulls_filler_pos+1] == -1))
		fulls_filler_pos += 2;

	while (partials_filler_pos < 2 * Q3C_NPARTIALS &&
		   !(partials_padded[partials_filler_pos] == 1 && partials_padded[partials_filler_pos+1] == -1))
		partials_filler_pos += 2;

	// Allocate memory for the arrays; allocate at least one byte so that a
	// zero-range result never hands NumPy a NULL data pointer.
	fulls    = malloc(fulls_filler_pos    > 0 ? fulls_filler_pos    * sizeof(q3c_ipix_t) : 1);
	partials = malloc(partials_filler_pos > 0 ? partials_filler_pos * sizeof(q3c_ipix_t) : 1);
	if (fulls == NULL || partials == NULL) {
		free(fulls);
		free(partials);
		return PyErr_NoMemory();
	}

	// copy the data since the padded arrays are on the caller's stack and will be the wrong length
	if (fulls_filler_pos > 0)
		memcpy(fulls, fulls_padded, fulls_filler_pos * sizeof(q3c_ipix_t));
	if (partials_filler_pos > 0)
		memcpy(partials, partials_padded, partials_filler_pos * sizeof(q3c_ipix_t));

	// Since the values being returned are ipix ranges, make it a 2D array.
	// data types: https://numpy.org/doc/stable/reference/c-api/dtype.html
	dims[0] = fulls_filler_pos/2;
	dims[1] = 2;
	np_fulls = (PyArrayObject *)PyArray_SimpleNewFromData(2, dims, NPY_INT64, fulls);
	if (np_fulls == NULL) {
		free(fulls);
		free(partials);
		return NULL;
	}
	PyArray_ENABLEFLAGS(np_fulls, NPY_ARRAY_OWNDATA); // set ownership of data: free when ndarray is garbage collected

	dims[0] = partials_filler_pos/2;
	np_partials = (PyArrayObject *)PyArray_SimpleNewFromData(2, dims, NPY_INT64, partials);
	if (np_partials == NULL) {
		Py_DECREF(np_fulls); // frees 'fulls' via OWNDATA
		free(partials);
		return NULL;
	}
	PyArray_ENABLEFLAGS(np_partials, NPY_ARRAY_OWNDATA);

	// Return as tuple. "N" steals the references; "O" would leak one reference
	// to each array per call (measured at ~890 bytes/call).
	return Py_BuildValue("NN", np_fulls, np_partials);
}

// Acquire read-only double buffers over a polygon's ra,dec vertex arrays and
// validate them. On success returns 0 with both views held (the caller must
// release them) and *n set to the vertex count; on failure returns -1 with an
// exception set and no views held.
//
// 'unwrapped_ra' is filled with the vertex ra values reduced to [0,360). Both
// polygon entry points reach q3c_get_facenum (through q3c_get_facenum_poly, which
// reads the first vertex), which requires that range; see the audit at the top of
// this file. A polygon written with negative ra - the natural way to describe one
// straddling ra=0 - otherwise picked the wrong cube face, and was then either
// projected wrongly or rejected as "too large". The values are copied rather than
// reduced in place: the buffer belongs to the caller's NumPy array.
// It must have room for Q3C_MAX_N_POLY_VERTEX values.
static int
polygon_vertex_views(PyObject *poly_ra, PyObject *poly_dec, Py_buffer *ra_view, Py_buffer *dec_view,
					 int *n, q3c_coord_t *unwrapped_ra)
{
	if (PyObject_GetBuffer(poly_ra, ra_view, PyBUF_FULL_RO) == -1)
		return -1;
	if (PyObject_GetBuffer(poly_dec, dec_view, PyBUF_FULL_RO) == -1) {
		PyBuffer_Release(ra_view);
		return -1;
	}

	#define POLY_VIEWS_FAIL() { PyBuffer_Release(ra_view); PyBuffer_Release(dec_view); return -1; }

	if (strcmp(ra_view->format, "d") != 0 || strcmp(dec_view->format, "d") != 0) {
		PyErr_SetString(PyExc_TypeError, "'poly_ra' and 'poly_dec' arrays must be of type np.double.");
		POLY_VIEWS_FAIL();
	}
	if (!PyBuffer_IsContiguous(ra_view, 'C') || !PyBuffer_IsContiguous(dec_view, 'C') ||
		ra_view->ndim != 1 || dec_view->ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "'poly_ra' and 'poly_dec' must be contiguous 1D arrays.");
		POLY_VIEWS_FAIL();
	}
	if (ra_view->shape[0] != dec_view->shape[0]) {
		PyErr_SetString(PyExc_ValueError, "'poly_ra' and 'poly_dec' arrays must be the same length.");
		POLY_VIEWS_FAIL();
	}
	*n = (int)ra_view->shape[0];
	if (*n < 3 || *n > Q3C_MAX_N_POLY_VERTEX) {
		PyErr_Format(PyExc_ValueError, "The polygon must have between 3 and %d vertices (was given %d).",
					 Q3C_MAX_N_POLY_VERTEX, *n);
		POLY_VIEWS_FAIL();
	}

	{
		const q3c_coord_t *ra_values = (const q3c_coord_t *)ra_view->buf;
		const q3c_coord_t *dec_values = (const q3c_coord_t *)dec_view->buf;
		int i;
		for (i = 0; i < *n; i++) {
			// a non-finite vertex reaches undefined float->int conversions in the
			// face selection below (upstream q3c issue #52)
			if (!isfinite(ra_values[i]) || !isfinite(dec_values[i])) {
				PyErr_SetString(PyExc_ValueError, "The polygon vertex values must all be finite.");
				POLY_VIEWS_FAIL();
			}
			unwrapped_ra[i] = UNWRAP_RA(ra_values[i]);
		}
	}

	#undef POLY_VIEWS_FAIL

	return 0;
}

static PyObject *
qlsc_q3c_radial_query(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	q3c_coord_t ra_cen, dec_cen, radius;

	// internal variables
	struct q3c_prm *hprm;
	q3c_ipix_t partials_padded[2 * Q3C_NPARTIALS];
	q3c_ipix_t fulls_padded[2 * Q3C_NFULLS];

	//static q3c_coord_t ra_cen_buf, dec_cen_buf, radius_buf;

	static char *kwlist[] = {"hprm", "ra", "dec", "radius", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "Oddd", // object + 2 doubles
									 kwlist,
									 &hprm_capsule, &ra_cen, &dec_cen, &radius))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	// Non-finite values pass the dec range check below (NaN > 90 is false) and
	// reach undefined float->int conversions (upstream q3c issue #52).
	if (!isfinite(ra_cen) || !isfinite(dec_cen) || !isfinite(radius))
	{
		PyErr_SetString(PyExc_ValueError, "The 'ra', 'dec', and 'radius' values must be finite.");
		return NULL;
	}

	ra_cen = UNWRAP_RA(ra_cen);
	if (q3c_fabs(dec_cen) > 90)
	{
		// the Python side validates this too; fail loudly if called directly
		// (PyErr_Format does not support float conversions, hence snprintf)
		char msg[128];
		snprintf(msg, sizeof msg, "The 'dec' value must be in the range [-90,90] (was given %.17g).", dec_cen);
		PyErr_SetString(PyExc_ValueError, msg);
		return NULL;
	}

	// generate the full, partials arrays
	// as of q3c v2.0.5, a nonzero return value means the ipix range arrays overflowed
	init_padded_ranges(fulls_padded, partials_padded);
	if (q3c_radial_query(hprm, ra_cen, dec_cen, radius, fulls_padded, partials_padded))
	{
		// %.17g preserves full double precision so the report is exactly reproducible
		char query[256];
		snprintf(query, sizeof query, "ra=%.17g, dec=%.17g, radius=%.17g "
				 "(degrees; the ra value is normalized to [0,360])", ra_cen, dec_cen, radius);
		raise_query_overflow("q3c_radial_query", query);
		return NULL;
	}

	return ipix_range_arrays(fulls_padded, partials_padded);
}

static PyObject *
qlsc_q3c_version(PyObject *module, PyObject *noargs)
{
	char version[32];
	q3c_get_version(version, sizeof version);
	version[sizeof version - 1] = '\0'; // q3c_get_version uses strncpy, which does not guarantee termination
	return PyUnicode_FromString(version);
}

static PyObject *
qlsc_q3c_check_sphere_point_in_poly(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	q3c_coord_t ra0, dec0;
	PyObject *np_poly_ra, *np_poly_dec;

	// internal variables
	struct q3c_prm *hprm;
	// per-face projection buffers: the main face plus up to four secondary faces.
	// NOTE: q3c proper (q3c.c) sizes these [3][...], which overflows for polygons
	// whose bounding box crosses onto three or four extra faces (multi_flag >= 3);
	// see the q3c_poly.c points[8] note in CLAUDE.md.
	q3c_coord_t xpj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t ypj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t poly_ra_unwrapped[Q3C_MAX_N_POLY_VERTEX]; // vertex ra reduced to [0,360)
	q3c_coord_t axpj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t aypj[5][Q3C_MAX_N_POLY_VERTEX];
	char faces[5];
	char multi_flag = 0;
	char too_large = 0;
	int n, result;
	Py_buffer ra_view, dec_view;

	static char *kwlist[] = {"hprm", "ra", "dec", "poly_ra", "poly_dec", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "OddOO", // object + 2 doubles + 2 objects
									 kwlist,
									 &hprm_capsule, &ra0, &dec0, &np_poly_ra, &np_poly_dec))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	// a non-finite point reaches an undefined float->int conversion inside
	// q3c_ang2ipix_xy (upstream q3c issue #52)
	if (!isfinite(ra0) || !isfinite(dec0)) {
		PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' values must be finite.");
		return NULL;
	}

	if (polygon_vertex_views(np_poly_ra, np_poly_dec, &ra_view, &dec_view, &n, poly_ra_unwrapped) != 0)
		return NULL;

	result = q3c_check_sphere_point_in_poly(hprm, n,
											poly_ra_unwrapped, (q3c_coord_t *)dec_view.buf,
											ra0, dec0, &too_large,
											0, // invocation = 0: always (re)project the polygon
											xpj, ypj, axpj, aypj, faces, &multi_flag);

	PyBuffer_Release(&ra_view);
	PyBuffer_Release(&dec_view);

	if (too_large) {
		PyErr_SetString(PyExc_ValueError, "The polygon is too large for the Q3C face projection; split it into smaller polygons.");
		return NULL;
	}

	return PyBool_FromLong(result);
}

static PyObject *
qlsc_q3c_check_sphere_points_in_poly(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	PyObject *np_ra, *np_dec;
	PyObject *np_poly_ra, *np_poly_dec;

	// internal variables
	struct q3c_prm *hprm;
	// per-face projection buffers; see the note in the scalar version above for why
	// there are five rows rather than q3c's three
	q3c_coord_t xpj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t ypj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t axpj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t aypj[5][Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t poly_ra_unwrapped[Q3C_MAX_N_POLY_VERTEX]; // vertex ra reduced to [0,360)
	char faces[5];
	char multi_flag = 0;
	char too_large = 0;
	int n;
	Py_buffer poly_ra_view, poly_dec_view;
	PyArrayObject *ra_array = NULL, *dec_array = NULL, *result_array = NULL;

	static char *kwlist[] = {"hprm", "ra", "dec", "poly_ra", "poly_dec", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "OOOOO", // 5 objects
									 kwlist,
									 &hprm_capsule, &np_ra, &np_dec, &np_poly_ra, &np_poly_dec))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	// PyArray_FROM_OTF deals with strides, non-contiguous input and integer types in
	// one call, as the ipix2* helpers above do
	ra_array = (PyArrayObject *)PyArray_FROM_OTF(np_ra, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (ra_array == NULL)
		return NULL;
	dec_array = (PyArrayObject *)PyArray_FROM_OTF(np_dec, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (dec_array == NULL) {
		Py_DECREF(ra_array);
		return NULL;
	}
	if (PyArray_NDIM(ra_array) != 1 || PyArray_NDIM(dec_array) != 1 ||
		PyArray_SIZE(ra_array) != PyArray_SIZE(dec_array)) {
		PyErr_SetString(PyExc_ValueError, "'ra' and 'dec' must be 1D arrays of the same length.");
		Py_DECREF(ra_array);
		Py_DECREF(dec_array);
		return NULL;
	}

	if (polygon_vertex_views(np_poly_ra, np_poly_dec, &poly_ra_view, &poly_dec_view, &n, poly_ra_unwrapped) != 0) {
		Py_DECREF(ra_array);
		Py_DECREF(dec_array);
		return NULL;
	}

	{
		const q3c_coord_t *ra_data = (const q3c_coord_t *)PyArray_DATA(ra_array);
		const q3c_coord_t *dec_data = (const q3c_coord_t *)PyArray_DATA(dec_array);
		q3c_coord_t *poly_dec = (q3c_coord_t *)poly_dec_view.buf;
		npy_intp n_points = PyArray_SIZE(ra_array);
		npy_intp i;
		npy_bool *result_data;
		npy_intp dims[1];

		for (i = 0; i < n_points; i++) {
			// a non-finite point reaches an undefined float->int conversion inside
			// q3c_ang2ipix_xy (upstream q3c issue #52)
			if (!isfinite(ra_data[i]) || !isfinite(dec_data[i])) {
				PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' arrays must contain only finite values.");
				goto points_in_poly_fail;
			}
		}

		dims[0] = n_points;
		result_array = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_BOOL);
		if (result_array == NULL)
			goto points_in_poly_fail;
		result_data = (npy_bool *)PyArray_DATA(result_array);

		// The polygon is projected onto the cube faces ONCE, by the invocation=0 call;
		// every later call passes invocation=1 and reuses the projection held in
		// xpj/ypj/axpj/aypj/faces/multi_flag. That is the whole point of this function:
		// the scalar wrapper reprojects all n vertices for every point tested, which is
		// the dominant cost of filtering a large candidate set through a polygon.
		// An empty point array still runs the projection, using the polygon's own first
		// vertex as the probe point, so that a polygon too large to project is reported
		// rather than passing silently.
		{
			q3c_coord_t probe_ra = (n_points > 0) ? ra_data[0] : poly_ra_unwrapped[0];
			q3c_coord_t probe_dec = (n_points > 0) ? dec_data[0] : poly_dec[0];
			int first = q3c_check_sphere_point_in_poly(hprm, n, poly_ra_unwrapped, poly_dec,
													   probe_ra, probe_dec, &too_large,
													   0, // invocation = 0: project the polygon
													   xpj, ypj, axpj, aypj, faces, &multi_flag);
			if (too_large) {
				PyErr_SetString(PyExc_ValueError, "The polygon is too large for the Q3C face projection; "
												  "split it into smaller polygons.");
				goto points_in_poly_fail;
			}
			if (n_points > 0)
				result_data[0] = first ? NPY_TRUE : NPY_FALSE;
		}

		for (i = 1; i < n_points; i++) {
			int inside = q3c_check_sphere_point_in_poly(hprm, n, poly_ra_unwrapped, poly_dec,
														ra_data[i], dec_data[i], &too_large,
														1, // invocation = 1: reuse the projection
														xpj, ypj, axpj, aypj, faces, &multi_flag);
			result_data[i] = inside ? NPY_TRUE : NPY_FALSE;
		}
	}

	PyBuffer_Release(&poly_ra_view);
	PyBuffer_Release(&poly_dec_view);
	Py_DECREF(ra_array);
	Py_DECREF(dec_array);

	// "N" steals the reference; "O" would leak one per call
	return Py_BuildValue("N", result_array);

points_in_poly_fail:
	PyBuffer_Release(&poly_ra_view);
	PyBuffer_Release(&poly_dec_view);
	Py_DECREF(ra_array);
	Py_DECREF(dec_array);
	Py_XDECREF(result_array);
	return NULL;
}

static PyObject *
qlsc_q3c_poly_query(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	PyObject *np_poly_ra, *np_poly_dec;

	// internal variables
	struct q3c_prm *hprm;
	q3c_ipix_t partials_padded[2 * Q3C_NPARTIALS];
	q3c_ipix_t fulls_padded[2 * Q3C_NFULLS];
	// The per-face projection of the polygon. Unlike q3c_check_sphere_point_in_poly,
	// which keeps one projection per face at once, q3c_poly_query reprojects into
	// these same arrays as it walks the faces, so one row each is right; they only
	// have to hold one vertex per polygon vertex.
	q3c_coord_t x[Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t y[Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t ax[Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t ay[Q3C_MAX_N_POLY_VERTEX];
	q3c_coord_t poly_ra_unwrapped[Q3C_MAX_N_POLY_VERTEX]; // vertex ra reduced to [0,360)
	q3c_poly qp;
	char too_large = 0;
	int n, overflowed;
	Py_buffer ra_view, dec_view;

	static char *kwlist[] = {"hprm", "poly_ra", "poly_dec", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "OOO", // 3 objects
									 kwlist,
									 &hprm_capsule, &np_poly_ra, &np_poly_dec))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	if (polygon_vertex_views(np_poly_ra, np_poly_dec, &ra_view, &dec_view, &n, poly_ra_unwrapped) != 0)
		return NULL;

	// q3c_project_poly only reads ra,dec; the projections are written to x,y,ax,ay
	qp.n = n;
	qp.ra  = poly_ra_unwrapped;
	qp.dec = (q3c_coord_t *)dec_view.buf;
	qp.x = x;
	qp.y = y;
	qp.ax = ax;
	qp.ay = ay;

	// q3c_poly_query returns early without filling the output arrays when the
	// polygon is too large, so they are initialized here first
	init_padded_ranges(fulls_padded, partials_padded);
	overflowed = q3c_poly_query(hprm, &qp, fulls_padded, partials_padded, &too_large);

	PyBuffer_Release(&ra_view);
	PyBuffer_Release(&dec_view);

	if (overflowed) {
		char query[64];
		snprintf(query, sizeof query, "the given polygon (%d vertices)", n);
		raise_query_overflow("q3c_poly_query", query);
		return NULL;
	}
	// checked before the arrays are read: too_large means they were never filled
	if (too_large) {
		PyErr_SetString(PyExc_ValueError, "The polygon is too large for the Q3C face projection; "
										  "split it into smaller polygons.");
		return NULL;
	}

	return ipix_range_arrays(fulls_padded, partials_padded);
}

// Validate the parameters shared by the ellipse functions. Returns 0 on success,
// -1 with an exception set otherwise. Note that a non-finite or out of range
// axis_ratio would otherwise reach q3c as a NaN eccentricity (sqrt of a negative),
// and non-finite coordinates reach undefined float->int conversions inside q3c
// (upstream issue #52).
static int
validate_ellipse_parameters(q3c_coord_t dec_cen, q3c_coord_t semi_major_axis,
							q3c_coord_t axis_ratio, q3c_coord_t position_angle)
{
	char msg[160];

	if (!isfinite(semi_major_axis) || !isfinite(axis_ratio) || !isfinite(position_angle)) {
		PyErr_SetString(PyExc_ValueError, "The 'semi_major_axis', 'axis_ratio' and 'position_angle' "
										  "values must be finite.");
		return -1;
	}
	if (q3c_fabs(dec_cen) > 90) {
		// (PyErr_Format does not support float conversions, hence snprintf)
		snprintf(msg, sizeof msg, "The 'dec' value must be in the range [-90,90] (was given %.17g).", dec_cen);
		PyErr_SetString(PyExc_ValueError, msg);
		return -1;
	}
	if (semi_major_axis < 0) {
		snprintf(msg, sizeof msg, "The 'semi_major_axis' value must not be negative (was given %.17g).",
				 semi_major_axis);
		PyErr_SetString(PyExc_ValueError, msg);
		return -1;
	}
	if (!(axis_ratio > 0 && axis_ratio <= 1)) {
		snprintf(msg, sizeof msg, "The 'axis_ratio' value must be in the range (0,1] (was given %.17g).",
				 axis_ratio);
		PyErr_SetString(PyExc_ValueError, msg);
		return -1;
	}
	// Q3C's ellipse is a conic in the gnomonic (tangent plane) projection about its
	// center, with semi-axes tan(semi_major_axis) and tan(semi_major_axis)*axis_ratio.
	// Only the near hemisphere projects, so the shape is meaningful strictly inside 90
	// degrees. Past it q3c_in_ellipse answers about a different region: it rejects every
	// point more than 90 degrees from the center whatever the semi-major axis, and its
	// sin^2 term stops being monotonic, so a 150 degree "ellipse" is smaller than a 60
	// degree one.
	//
	// This is checked for BOTH exported ellipse functions, the predicate and the
	// candidate range query, though only the predicate can return a wrong answer today:
	// an oversized range query is saved by upstream's Q3C_MAXRAD fallback, which returns
	// the whole sphere and so is a safe superset. But that safety is a fact about
	// q3cube.c rather than about this wrapper, and an upstream sync that dropped the
	// fallback would silently start under-covering a region that was not the one asked
	// for. Refusing here costs nothing - the Python layer answers a circle (axis_ratio
	// == 1) of any radius with a cone search instead, and never asks for an ellipse this
	// large. The fallback itself is pinned by source/tests/test_q3c_maxrad_fallback.py so
	// a sync cannot remove it unnoticed either.
	if (!(semi_major_axis < ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS)) {
		char limit_msg[320];
		snprintf(limit_msg, sizeof limit_msg,
				 "The 'semi_major_axis' value must be less than %g degrees (was given %.17g): Q3C's ellipse is "
				 "a conic in the tangent plane at its center, which only reaches the near hemisphere. Use "
				 "axis_ratio=1 (a circle) or radial_query for a larger region.",
				 (double)ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS, semi_major_axis);
		PyErr_SetString(PyExc_ValueError, limit_msg);
		return -1;
	}
	return 0;
}

// q3c takes the eccentricity, while the value everything else here (and q3c's own
// SQL interface) speaks in is the axis ratio b/a. This is the conversion q3c.c makes.
static q3c_coord_t
ellipse_eccentricity(q3c_coord_t axis_ratio)
{
	return q3c_sqrt(1 - axis_ratio * axis_ratio);
}

static PyObject *
qlsc_q3c_ellipse_query(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	q3c_coord_t ra_cen, dec_cen, semi_major_axis, axis_ratio, position_angle;

	// internal variables
	struct q3c_prm *hprm;
	q3c_ipix_t partials_padded[2 * Q3C_NPARTIALS];
	q3c_ipix_t fulls_padded[2 * Q3C_NFULLS];

	static char *kwlist[] = {"hprm", "ra", "dec", "semi_major_axis", "axis_ratio", "position_angle", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "Oddddd", // object + 5 doubles
									 kwlist,
									 &hprm_capsule, &ra_cen, &dec_cen,
									 &semi_major_axis, &axis_ratio, &position_angle))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

	if (!isfinite(ra_cen) || !isfinite(dec_cen)) {
		PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' values must be finite.");
		return NULL;
	}
	if (validate_ellipse_parameters(dec_cen, semi_major_axis, axis_ratio, position_angle) != 0)
		return NULL;

	ra_cen = UNWRAP_RA(ra_cen);

	// as of q3c v2.0.5, a nonzero return value means the ipix range arrays overflowed
	init_padded_ranges(fulls_padded, partials_padded);
	if (q3c_ellipse_query(hprm, ra_cen, dec_cen, semi_major_axis,
						  ellipse_eccentricity(axis_ratio), position_angle,
						  fulls_padded, partials_padded))
	{
		// %.17g preserves full double precision so the report is exactly reproducible
		char query[384];
		snprintf(query, sizeof query,
				 "ra=%.17g, dec=%.17g, semi_major_axis=%.17g, axis_ratio=%.17g, position_angle=%.17g "
				 "(degrees; the ra value is normalized to [0,360])",
				 ra_cen, dec_cen, semi_major_axis, axis_ratio, position_angle);
		raise_query_overflow("q3c_ellipse_query", query);
		return NULL;
	}

	return ipix_range_arrays(fulls_padded, partials_padded);
}

static PyObject *
qlsc_q3c_in_ellipse(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t ra, dec;                 // the point being tested
	q3c_coord_t ra_cen, dec_cen;         // the center of the ellipse
	q3c_coord_t semi_major_axis, axis_ratio, position_angle;

	// internal variables
	char result;

	static char *kwlist[] = {"ra", "dec", "ra_cen", "dec_cen",
							 "semi_major_axis", "axis_ratio", "position_angle", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "ddddddd", // 7 doubles
									 kwlist,
									 &ra, &dec, &ra_cen, &dec_cen,
									 &semi_major_axis, &axis_ratio, &position_angle))
	{
		// unable to parse inputs -> raise exception
		return NULL;
	}

	if (!isfinite(ra) || !isfinite(dec) || !isfinite(ra_cen) || !isfinite(dec_cen)) {
		PyErr_SetString(PyExc_ValueError, "The 'ra' and 'dec' values must be finite.");
		return NULL;
	}
	if (validate_ellipse_parameters(dec_cen, semi_major_axis, axis_ratio, position_angle) != 0)
		return NULL;
	if (q3c_fabs(dec) > 90) {
		char msg[160];
		snprintf(msg, sizeof msg, "The 'dec' value must be in the range [-90,90] (was given %.17g).", dec);
		PyErr_SetString(PyExc_ValueError, msg);
		return NULL;
	}
	// argument order as in q3c.c's pgq3c_in_ellipse: center first, then the point
	result = q3c_in_ellipse(ra_cen, dec_cen, ra, dec, semi_major_axis,
							ellipse_eccentricity(axis_ratio), position_angle);

	return PyBool_FromLong(result);
}

static PyMethodDef qlsc_methods[] = { // METH_VARARGS _or_ METH_VARARGS | METH_KEYWORDS
	// Ref: available flags: https://docs.python.org/3/c-api/structures.html#c.PyMethodDef
	
	// This is an array of PyMethodDef structs, NULL-terminated
	// Format: {function name in Python, pointer to static C function, flags indicating inputs, docstring}
	
	// flags:
	// METH_VARARGS -> f(self, args)
	// METH_NOARGS -> no arguments
	// METH_O -> one argument
	// METH_VARARGS | METH_KEYWORDS -> f(self, args, kwargs )
	{"init_q3c", (PyCFunction)qlsc_init_q3c1, METH_VARARGS|METH_KEYWORDS, "Initialize prm, Q3C's main structure."},
	{"nside", (PyCFunction)qlsc_q3c_nside, METH_VARARGS, "Return the number of bins along the edge of each cube face."},
	{"version", (PyCFunction)qlsc_q3c_version, METH_NOARGS, "Return the version of the embedded Q3C code."},
	{"ang2ipix", (PyCFunction)qlsc_q3c_ang2ipix, METH_VARARGS|METH_KEYWORDS, "Convert ra,dec to ipix value."},
	{"ang2ipix_xy", (PyCFunction)qlsc_q3c_ang2ipix_xy, METH_VARARGS|METH_KEYWORDS, "Convert ra,dec to ipix value, also returning (x,y) and the face number in a dictionary."},
	{"ipix2ang", (PyCFunction)qlsc_q3c_ipix2ang, METH_VARARGS|METH_KEYWORDS, "Convert an ipix value (or an array of them) to a ra,dec tuple (or a tuple of ra,dec arrays)."},
	{"ipix2xy", (PyCFunction)qlsc_q3c_ipix2xy, METH_VARARGS|METH_KEYWORDS, "Convert an ipix value (or an array of them) to the (x,y) coordinate on the square face with the face number as a tuple: (facenum,x,y); array input returns a tuple of three arrays."},
	{"facenum", (PyCFunction)qlsc_q3c_facenum, METH_VARARGS|METH_KEYWORDS, "Return the cube face number for the given coordinates."},
	{"pixarea", (PyCFunction)qlsc_q3c_pixarea, METH_VARARGS|METH_KEYWORDS, "Return the area of a given Q3C pixel for a given ipix and depth."},
	{"distance", (PyCFunction)qlsc_q3c_dist, METH_VARARGS|METH_KEYWORDS, "Calculates the angular distance between two points on a sphere, in degrees."},
	{"sindist", (PyCFunction)qlsc_q3c_sindist, METH_VARARGS|METH_KEYWORDS, "Calculates the haversine sin^2(d/2) of the angular distance d between two points on a sphere - NOT sin(d); see qlsc.sindist."},
	{"xy2ang", (PyCFunction)qlsc_q3c_xy2ang, METH_VARARGS|METH_KEYWORDS, "Convert an x,y coordinate pair on the given face number to (ra,dec)."},
	{"xy2facenum", (PyCFunction)qlsc_q3c_xy2facenum, METH_VARARGS|METH_KEYWORDS, "Convert an x,y coordinate pair on the given face number to the corresponding cube face number."},
	{"check_sphere_point_in_poly", (PyCFunction)qlsc_q3c_check_sphere_point_in_poly, METH_VARARGS|METH_KEYWORDS, "Test whether an ra,dec point lies inside a spherical polygon."},
	{"check_sphere_points_in_poly", (PyCFunction)qlsc_q3c_check_sphere_points_in_poly, METH_VARARGS|METH_KEYWORDS, "Test an array of ra,dec points against a spherical polygon, projecting the polygon once; returns a boolean array."},
	{"in_ellipse", (PyCFunction)qlsc_q3c_in_ellipse, METH_VARARGS|METH_KEYWORDS, "Test whether an ra,dec point lies inside an ellipse on the sphere."},
	{"radial_query", (PyCFunction)qlsc_q3c_radial_query, METH_VARARGS|METH_KEYWORDS, "Return the (fulls,partials) ipix ranges covering a cone."},
	{"poly_query", (PyCFunction)qlsc_q3c_poly_query, METH_VARARGS|METH_KEYWORDS, "Return the (fulls,partials) ipix ranges covering a spherical polygon."},
	{"ellipse_query", (PyCFunction)qlsc_q3c_ellipse_query, METH_VARARGS|METH_KEYWORDS, "Return the (fulls,partials) ipix ranges covering an ellipse on the sphere."},
	{NULL, NULL, 0, NULL}	// sentinel
};

// Ref: https://docs.python.org/3/c-api/module.html#c.PyModuleDef
static struct PyModuleDef qlsc_module_definition = {
	PyModuleDef_HEAD_INIT,
	"qlsc_c",					// name of module
	"A Python implementation of the Quadrilateralized Spherical Cube (QLSC).",	// docstring for the module
	-1,						// size of per-interpreter state of the module,
			                // or -1 if the module keeps state in global variables.
	qlsc_methods			// pointer to a table of module-level functions (PyMethodDef class)
//	NULL,					// m_slots
//	NULL,					// m_traverse
//	NULL,					// m_clear, a clear function to call during GC clearing of the module object
//	NULL					// m_free, function to call during deallocation of the module object
};

PyMODINIT_FUNC
PyInit_q3c(void) // must be named "PyInit_" + name of extension (without package name)
{
	import_array(); // initialize NumPy, see: https://numpy.org/doc/stable/reference/c-api/array.html#importing-the-api
	PyObject *module = PyModule_Create(&qlsc_module_definition);
	if (module == NULL)
		return NULL;

	// a dedicated exception type for q3c's internal range-capacity overflow, so
	// callers can catch it without matching on a generic RuntimeError
	QlscRangeOverflowError = PyErr_NewException("qlsc.q3c.RangeOverflowError", PyExc_RuntimeError, NULL);
	if (QlscRangeOverflowError == NULL || PyModule_AddObject(module, "RangeOverflowError", QlscRangeOverflowError) < 0) {
		Py_XDECREF(QlscRangeOverflowError);
		Py_DECREF(module);
		return NULL;
	}

	// Two limits the Python layer has to reason about, published rather than
	// duplicated as literals up there:
	//
	//   Q3C_MAXRAD - the size at and above which the vendored query functions stop
	//     doing geometry and return the whole sphere. Growing a search past it cannot
	//     widen the coverage any further, so it is the ceiling for the widening loop.
	//   ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS - the largest elliptical region this wrapper
	//     will accept; see validate_ellipse_parameters.
	//
	// Both are added as floats rather than with PyModule_AddIntConstant: Q3C_MAXRAD is
	// an integer literal today, but a future upstream release could make it fractional
	// and an int conversion would silently truncate it.
	{
		PyObject *maxrad = PyFloat_FromDouble((double)Q3C_MAXRAD);
		if (maxrad == NULL || PyModule_AddObject(module, "Q3C_MAXRAD", maxrad) < 0) {
			Py_XDECREF(maxrad);
			Py_DECREF(module);
			return NULL;
		}
		PyObject *max_semi_major = PyFloat_FromDouble((double)ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS);
		if (max_semi_major == NULL ||
			PyModule_AddObject(module, "ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS", max_semi_major) < 0) {
			Py_XDECREF(max_semi_major);
			Py_DECREF(module);
			return NULL;
		}
	}

	return module;
}


