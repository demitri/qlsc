#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/ndarrayobject.h> // needed to work with NumPy arrays

#include "q3c/common.h"
//#include "my_bits.h"

// character parameter values
// ref: https://docs.python.org/3/c-api/arg.html?highlight=pyarg_parsetupleandkeywords#building-values

#define Q3C_STRUCT_POINTER_BUFFER "Q3C_prm_struct_pointer"

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
				//PySys_WriteStdout("ra,dec = (%f,%f)\n", ra, dec);
				q3c_ang2ipix(hprm, ra, dec, &ipix_array[i]);
			}
			
		} else {
			// standard C array, contiguous data
			
			ipix_array = malloc(n * sizeof(q3c_ipix_t)); // allocate the ipix values array
			
			for (int i=0; i < n; i++) {
				q3c_ang2ipix(hprm, ra_buffer[i], dec_buffer[i], &ipix_array[i]);
				//PySys_WriteStdout("C ra,dec = (%f,%f) -> %lld\n", ra_buffer[i], dec_buffer[i], ipix_array[i]);
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

static PyObject *
qlsc_q3c_ipix2ang(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
{
	// external parameters
	q3c_ipix_t ipix; // int64_t
	struct q3c_prm *hprm;
	
	// internal variables
	q3c_coord_t ra; // double
	q3c_coord_t dec;
	PyObject *hprm_capsule;

	static char *kwlist[] = {"hprm", "ipix", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "OL",			// object + longlong
    								 kwlist,
    								 &hprm_capsule, &ipix))
	{
		// unable to parse inputs -> raise exception
        return NULL;
	}
		
	//PyObject *hprm_attr_name = PyUnicode_FromString("_hprm");
	//PySys_WriteStdout("has attr: %d\n", PyObject_HasAttr(self, hprm_attr_name));

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception

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
	q3c_ipix_t *xbits1;
	q3c_ipix_t *ybits1;
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
	q3c_ipix_t ipix;
	PyObject *hprm_capsule;
	
	// internal variables
	struct q3c_prm *hprm;
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
	
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "OL",						// object + longlong
    								 kwlist,					//
    								 &hprm_capsule, &ipix))
	{
		// unable to parse inputs -> raise exception
        return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	if (hprm == NULL)
		return NULL; // not a valid hprm capsule; PyCapsule_GetPointer has set the exception
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

// Raise a RangeOverflowError for a q3c_radial_query ipix range overflow,
// reporting the query values as used (i.e. after ra normalization).
static void raise_radial_query_overflow(q3c_coord_t ra, q3c_coord_t dec, q3c_coord_t radius)
{
	char msg[512];
	// %.17g preserves full double precision so the report is exactly reproducible
	snprintf(msg, sizeof msg,
			 "q3c_radial_query: too many ipix ranges for the query ra=%.17g, dec=%.17g, radius=%.17g "
			 "(degrees; the ra value is normalized to [0,360]). This is an internal error in the embedded "
			 "Q3C code (version " Q3C_VERSION "); please report it with these values at "
			 "https://github.com/segasai/q3c/issues .",
			 ra, dec, radius);
	PyErr_SetString(QlscRangeOverflowError, msg);
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
	q3c_ipix_t *partials, *fulls;
	int fulls_filler_pos, partials_filler_pos, max_idx;
	npy_intp dims[2]; // aka const long *
	int fulls_length = 2 * Q3C_NFULLS;
	int partials_length = 2 * Q3C_NPARTIALS;
	
	// returned values
	PyArrayObject *np_fulls;
	PyArrayObject *np_partials;

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
	if (q3c_radial_query(hprm, ra_cen, dec_cen, radius, fulls_padded, partials_padded))
	{
		raise_radial_query_overflow(ra_cen, dec_cen, radius);
		return NULL;
	}

	// q3c_radial_query() calls array_filler() which adds (1,-1)
	// pairs to fill out the arrays. We don't need them here, so we'll
	// find number of actual data points and dynamically create arrays of the needed size.
	// Only the exact (1,-1) pair is filler - a legitimate range can start or end at 1.
	fulls_filler_pos = 0;
	max_idx = fulls_length;
	while (fulls_filler_pos < max_idx && !(fulls_padded[fulls_filler_pos] == 1 && fulls_padded[fulls_filler_pos+1] == -1))
		fulls_filler_pos += 2;

	partials_filler_pos = 0;
	max_idx = partials_length;
	while (partials_filler_pos < max_idx && !(partials_padded[partials_filler_pos] == 1 && partials_padded[partials_filler_pos+1] == -1))
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

	// copy the data since the padded arrays are on the stack and will be the wrong length
	if (fulls_filler_pos > 0)
		memcpy(fulls, fulls_padded, fulls_filler_pos * sizeof(q3c_ipix_t));
	if (partials_filler_pos > 0)
		memcpy(partials, partials_padded, partials_filler_pos * sizeof(q3c_ipix_t));

	// Create NumPy arrays and return them.
	// Since the values being returned are ipix ranges, make it a 2D array.
	// data types: https://numpy.org/doc/stable/reference/c-api/dtype.html
	//
	dims[0] = fulls_filler_pos/2;
	dims[1] = 2;
	np_fulls = (PyArrayObject *)PyArray_SimpleNewFromData(2, 			// int nd - length of dims array
														  dims, 		// npy_intp* dims (for 1D, could just be an int)
														  NPY_INT64, 	// int typenum,
														  fulls); 		// void* data
	if (np_fulls == NULL) {
		free(fulls);
		free(partials);
		return NULL;
	}
	PyArray_ENABLEFLAGS(np_fulls, NPY_ARRAY_OWNDATA); // set ownership of data: free when ndarray is garbage collected

	dims[0] = partials_filler_pos/2;
	np_partials = (PyArrayObject *)PyArray_SimpleNewFromData(2, 			// int nd - length of dims array
															 dims, 			// npy_intp* dims (for 1D, could just be an int)
															 NPY_INT64, 	// int typenum,
															 partials); 	// void* data
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

	if (PyObject_GetBuffer(np_poly_ra, &ra_view, PyBUF_FULL_RO) == -1)
		return NULL;
	if (PyObject_GetBuffer(np_poly_dec, &dec_view, PyBUF_FULL_RO) == -1) {
		PyBuffer_Release(&ra_view);
		return NULL;
	}

	// the checks below jump here on error
	#define POLY_VIEWS_FAIL() { PyBuffer_Release(&ra_view); PyBuffer_Release(&dec_view); return NULL; }

	if (strcmp(ra_view.format, "d") != 0 || strcmp(dec_view.format, "d") != 0) {
		PyErr_SetString(PyExc_TypeError, "'poly_ra' and 'poly_dec' arrays must be of type np.double.");
		POLY_VIEWS_FAIL();
	}
	if (!PyBuffer_IsContiguous(&ra_view, 'C') || !PyBuffer_IsContiguous(&dec_view, 'C') ||
		ra_view.ndim != 1 || dec_view.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "'poly_ra' and 'poly_dec' must be contiguous 1D arrays.");
		POLY_VIEWS_FAIL();
	}
	if (ra_view.shape[0] != dec_view.shape[0]) {
		PyErr_SetString(PyExc_ValueError, "'poly_ra' and 'poly_dec' arrays must be the same length.");
		POLY_VIEWS_FAIL();
	}
	n = (int)ra_view.shape[0];
	if (n < 3 || n > Q3C_MAX_N_POLY_VERTEX) {
		PyErr_Format(PyExc_ValueError, "The polygon must have between 3 and %d vertices (was given %d).",
					 Q3C_MAX_N_POLY_VERTEX, n);
		POLY_VIEWS_FAIL();
	}

	#undef POLY_VIEWS_FAIL

	result = q3c_check_sphere_point_in_poly(hprm, n,
											(q3c_coord_t *)ra_view.buf, (q3c_coord_t *)dec_view.buf,
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
	{"ipix2ang", (PyCFunction)qlsc_q3c_ipix2ang, METH_VARARGS|METH_KEYWORDS, "Convert an ipix value to ra,dec tuple."},
	{"ipix2xy", (PyCFunction)qlsc_q3c_ipix2xy, METH_VARARGS|METH_KEYWORDS, "Convert an ipix value to the (x,y) coordinate on the square face with the face number as a tuple: (facenum,x,y)."},
	{"facenum", (PyCFunction)qlsc_q3c_facenum, METH_VARARGS|METH_KEYWORDS, "Return the cube face number for the given coordinates."},
	{"pixarea", (PyCFunction)qlsc_q3c_pixarea, METH_VARARGS|METH_KEYWORDS, "Return the area of a given Q3C pixel for a given ipix and depth."},
	{"distance", (PyCFunction)qlsc_q3c_dist, METH_VARARGS|METH_KEYWORDS, "Calculates angular distance between two points on a sphere."},
	{"sindist", (PyCFunction)qlsc_q3c_sindist, METH_VARARGS|METH_KEYWORDS, "Calculates the sine of the angular distance between two points on a sphere."},
	{"xy2ang", (PyCFunction)qlsc_q3c_xy2ang, METH_VARARGS|METH_KEYWORDS, "Convert an x,y coordinate pair on the given face number to (ra,dec)."},
	{"xy2facenum", (PyCFunction)qlsc_q3c_xy2facenum, METH_VARARGS|METH_KEYWORDS, "Convert an x,y coordinate pair on the given face number to the corresponding cube face number."},
	{"check_sphere_point_in_poly", (PyCFunction)qlsc_q3c_check_sphere_point_in_poly, METH_VARARGS|METH_KEYWORDS, "Test whether an ra,dec point lies inside a spherical polygon."},
	{"radial_query", (PyCFunction)qlsc_q3c_radial_query, METH_VARARGS|METH_KEYWORDS, ""},
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

	return module;
}


