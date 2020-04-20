#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "common.h"
//#include "my_bits.h"

// character parameter values
// ref: https://docs.python.org/3/c-api/arg.html?highlight=pyarg_parsetupleandkeywords#building-values

#define Q3C_STRUCT_POINTER_BUFFER "Q3C_prm_struct_pointer"

// ---------------------------------
// prm (main Q3C structure) routines
// ---------------------------------
// TODO: submit this as a pull request to include in q3cube.c ?
void free_q3c1(struct q3c_prm *hprm)
{
	free(hprm->xbits);
	free(hprm->ybits);
	free(hprm->xbits1);
	free(hprm->ybits1);
	free(hprm);
}

// Function that Python will call when the prm (wrapped in a PyObject) is freed.
static void del_prm(PyObject *obj)
{
	// get C pointer out of object
	struct q3c_prm *hprm = PyCapsule_GetPointer(obj, Q3C_STRUCT_POINTER_BUFFER);
	free_q3c1(hprm);
}

static PyObject *
pyq3c_init_q3c1(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *returnValue = NULL;
	
	// dynamically allocate memory for this as it will be kept around
	struct q3c_prm *hprm = malloc(sizeof (struct q3c_prm));
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
pyq3c_q3c_nside(PyObject *self, PyObject *arg)
{
	PyObject *return_value = NULL;

	struct q3c_prm *hprm;
	
	// no parsing needed; one value is expected
	//hprm = (struct q3c_prm*)PyCapsule_GetPointer(arg, Q3C_STRUCT_POINTER_BUFFER);
	
	PyObject *capsule;
	PyArg_ParseTuple(arg, "O", &capsule);
	hprm = PyCapsule_GetPointer(capsule, Q3C_STRUCT_POINTER_BUFFER);
	
	return_value = PyLong_FromLongLong(hprm->nside);
	
//	goto finally;
//except:
//	return_value = NULL;
//finally:
	return return_value;
}

static PyObject *
pyq3c_q3c_ang2ipix(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
//pyq3c_q3c_ang2ipix(PyObject *module, PyObject *args) // -> cast as PyCFunction
{
	// external parameters
	q3c_coord_t ra;
	q3c_coord_t dec;
	PyObject *hprm_capsule;
	
	// internal variables
	struct q3c_prm *hprm;
	q3c_ipix_t ipix = 0;
	static int invocation;
	static q3c_coord_t ra_buf, dec_buf;
	static q3c_ipix_t ipix_buf;

	static char *kwlist[] = {"hprm", "ra", "dec", NULL};

	//PySys_WriteStdout("about to parse input\n");
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
    								 "Odd",	// object + 2 doubles
    								 kwlist,
    								 &hprm_capsule, &ra, &dec))
	{
		// unable to parse inputs -> raise exception
		PySys_WriteStdout("unable to parse input, returning NULL\n");
        return NULL;
	}

	//PySys_WriteStdout("ra,dec = %.6f, %.6f\n", ra, dec);
	
	// TODO: look for examples where attribute names are pre-cached.
	//PyObject *hprm_attr_name = PyUnicode_FromString("_hprm");
	//PySys_WriteStdout("has attr: %d\n", PyObject_HasAttr(self, hprm_attr_name));

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	
	if (invocation==0)
	{
		;
	}
	else
	{
		if ((ra == ra_buf) && (dec == dec_buf))
		{
			PyLong_FromLongLong(ipix_buf);
		}
	}
	if ((!isfinite(ra)) || (!isfinite(dec)))
	{
		Py_RETURN_NONE;
	}
		
	q3c_ang2ipix(hprm, ra, dec, &ipix);

	ra_buf = ra;
	dec_buf = dec;
	ipix_buf = ipix;
	invocation=1;

	//PySys_WriteStdout("ipix = %" PRId64 "\n", ipix);

    return PyLong_FromLongLong(ipix);
}

static PyObject *
pyq3c_q3c_ang2ipix_xy(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
{
	// external parameters
	q3c_coord_t ra;
	q3c_coord_t dec;
	PyObject *hprm_capsule;

	//internal variables
	struct q3c_prm *hprm;
	static int invocation;
	static q3c_coord_t ra_buf, dec_buf;
	static q3c_ipix_t ipix_buf;

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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
        return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);

	if (invocation==0)
	{
		;
	}
	else
	{
		if ((ra == ra_buf) && (dec == dec_buf))
		{
			PyLong_FromLongLong(ipix_buf);
		}
	}
	if ((!isfinite(ra)) || (!isfinite(dec)))
	{
		Py_RETURN_NONE;
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
pyq3c_q3c_ipix2ang(PyObject *module, PyObject *args, PyObject *kwargs) // -> cast as PyCFunctionWithKeywords
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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
        return NULL;
	}
		
	//PyObject *hprm_attr_name = PyUnicode_FromString("_hprm");
	//PySys_WriteStdout("has attr: %d\n", PyObject_HasAttr(self, hprm_attr_name));

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);

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
pyq3c_q3c_ipix2xy(PyObject *module, PyObject *args, PyObject *kwargs)
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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
        return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	//face_num = ipix / (hprm->nside)^2;
	face_num = (int)trunc(ipix/pow(hprm->nside,2)); // take care of types here...

//	PySys_WriteStdout("x = %f\n", (double) (ipix / ((hprm->nside)^2)) );
//	PySys_WriteStdout("ipix = %"PRId64"\n", ipix);
//	PySys_WriteStdout("face_num = %d\n", face_num);
//	PySys_WriteStdout("hprm->nside = %f\n", pow(hprm->nside,2));
	
	ipix_to_xy(hprm, ipix, &x, &y);
	
	return Py_BuildValue("ldd", face_num, x, y);
}

static PyObject *
pyq3c_q3c_pixarea(PyObject *module, PyObject *args, PyObject *kwargs)
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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
        return NULL;
	}
	
	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	
	area = q3c_pixarea(hprm, ipix, depth);

	// debug
	//PySys_WriteStdout("ipix = %"PRId64"\n", ipix);
	//PySys_WriteStdout("depth = %d\n", depth);
	//PySys_WriteStdout("area = %.10e\n", area);
	
	return PyFloat_FromDouble(area);
}

static PyObject *
pyq3c_q3c_facenum(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	q3c_coord_t ra, dec;
	PyObject *hprm_capsule;
	
	// internal variables
	struct q3c_prm *hprm;
	char facenum;
	
	static char *kwlist[] = {"hprm", "ra", "dec", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									  "Odd",			// object + 2 doubles
									  kwlist,
									  &hprm_capsule, &ra, &dec))
	 {
		 // unable to parse inputs -> raise exception
		 PySys_WriteStdout("unable to parse input, returning NULL\n");
		 return NULL;
	 }
	
	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);
	
	facenum = q3c_get_facenum(ra, dec);
	
	return PyLong_FromLong((long)facenum);
}

static PyObject *
pyq3c_q3c_dist(PyObject *module, PyObject *args, PyObject *kwargs)
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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
		return NULL;
	}

	distance = q3c_dist( ra1, dec1, ra2, dec2 );
	
	return PyFloat_FromDouble(distance);
}

static PyObject *
pyq3c_q3c_sindist(PyObject *module, PyObject *args, PyObject *kwargs)
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
		// unable to parse inputs -> raise exception
		PySys_WriteStdout("unable to parse input, returning NULL\n");
		return NULL;
	}

	sine_distance = q3c_sindist( ra1, dec1, ra2, dec2 );
	
	return PyFloat_FromDouble(sine_distance);
}

static PyObject *
pyq3c_q3c_xy2ang(PyObject *module, PyObject *args, PyObject *kwargs)
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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
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

/* I don't really know what this does.
static PyObject *
pyq3c_q3c_xy2facenum(PyObject *module, PyObject *args, PyObject *kwargs)
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
		PySys_WriteStdout("unable to parse input, returning NULL\n");
		return NULL;
	}
	
	char result = q3c_xy2facenum(x, y, facenum);
	PySys_WriteStdout("(%.4f, %.4f) | %d -> %hhd\n", x, y, facenum, result);
		
	return PyLong_FromLong( (long) result );
}
*/

static PyObject *
pyq3c_q3c_radial_query_it(PyObject *module, PyObject *args, PyObject *kwargs)
{
	// external parameters
	PyObject *hprm_capsule;
	q3c_coord_t ra_cen, dec_cen, radius;
	int iteration;
	int full_flag; //  1 = full, 0 = partial
	
	// internal variables
	struct q3c_prm *hprm;
	static int invocation = 0;
	
	static q3c_coord_t ra_cen_buf, dec_cen_buf, radius_buf;
	
	static q3c_ipix_t partials[2 * Q3C_NPARTIALS];
	static q3c_ipix_t fulls[2 * Q3C_NFULLS];

	static char *kwlist[] = {"hprm", "ra", "dec", "radius", "iteration", "full_flag", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs,
									 "Odddii", // object + 2 doubles + 2 ints
									 kwlist,
									 &hprm_capsule, &ra_cen, &dec_cen, &radius, &iteration, &full_flag))
	{
		// unable to parse inputs -> raise exception
		PySys_WriteStdout("unable to parse input, returning NULL\n");
		return NULL;
	}

	hprm = (struct q3c_prm*)PyCapsule_GetPointer(hprm_capsule, Q3C_STRUCT_POINTER_BUFFER);

	ra_cen = UNWRAP_RA(ra_cen);
	if (q3c_fabs(dec_cen) > 90)
	{
		PySys_WriteStdout("'dec' value out of range - todo: raise exception'\n");
	}
	
	if (invocation == 0) {
		;
	} else {
		if ((ra_cen == ra_cen_buf) && (dec_cen == dec_cen_buf) &&
			radius == radius_buf)
		{
			if (full_flag)
				return PyLong_FromLongLong(fulls[iteration]);
			else
				return PyLong_FromLongLong(partials[iteration]);
		}
	}
	
	q3c_radial_query(hprm, ra_cen, dec_cen, radius, fulls, partials);
	
	// cache values
	ra_cen_buf = ra_cen;
	dec_cen_buf = dec_cen;
	radius_buf = radius;
	invocation = 1;
	
	return full_flag ? PyLong_FromLongLong(fulls[iteration]) : PyLong_FromLongLong(partials[iteration]);
}

static PyMethodDef pyq3c_methods[] = { // METH_VARARGS _or_ METH_VARARGS | METH_KEYWORDS
	// Ref: available flags: https://docs.python.org/3/c-api/structures.html#c.PyMethodDef
	
	// This is an array of PyMethodDef structs, NULL-terminated
	// Format: {function name in Python, pointer to static C function, flags indicating inputs, docstring}
	
	// flags:
	// METH_VARARGS -> f(self, args)
	// METH_NOARGS -> no arguments
	// METH_O -> one argument
	// METH_VARARGS | METH_KEYWORDS -> f(self, args, kwargs )
	{"init_q3c", (PyCFunction)pyq3c_init_q3c1, METH_VARARGS|METH_KEYWORDS, "Initialize prm, Q3C's main structure."},
	{"nside", (PyCFunction)pyq3c_q3c_nside, METH_VARARGS, "Return the number of bins along the edge of each cube face."},
	{"ang2ipix", (PyCFunction)pyq3c_q3c_ang2ipix, METH_VARARGS|METH_KEYWORDS, "Convert ra,dec to ipix value."},
	{"ang2ipix_xy", (PyCFunction)pyq3c_q3c_ang2ipix_xy, METH_VARARGS|METH_KEYWORDS, "Convert ra,dec to ipix value, also returning (x,y) and the face number in a dictionary."},
	{"ipix2ang", (PyCFunction)pyq3c_q3c_ipix2ang, METH_VARARGS|METH_KEYWORDS, "Convert an ipix value to ra,dec tuple."},
	{"ipix2xy", (PyCFunction)pyq3c_q3c_ipix2xy, METH_VARARGS|METH_KEYWORDS, "Convert an ipix value to the (x,y) coordinate on the square face with the face number as a tuple: (facenum,x,y)."},
	{"facenum", (PyCFunction)pyq3c_q3c_facenum, METH_VARARGS|METH_KEYWORDS, "Return the cube face number for the given coordinates."},
	{"pixarea", (PyCFunction)pyq3c_q3c_pixarea, METH_VARARGS|METH_KEYWORDS, "Return the area of a given Q3C pixel for a given ipix and depth."},
	{"distance", (PyCFunction)pyq3c_q3c_dist, METH_VARARGS|METH_KEYWORDS, "Calculates angular distance between two points on a sphere."},
	{"sindist", (PyCFunction)pyq3c_q3c_sindist, METH_VARARGS|METH_KEYWORDS, "Calculates the sine of the angular distance between two points on a sphere."},
	{"xy2ang", (PyCFunction)pyq3c_q3c_xy2ang, METH_VARARGS|METH_KEYWORDS, "Convert an x,y coordinate pair on the given face number to (ra,dec)."},
//	{"xy2facenum", (PyCFunction)pyq3c_q3c_xy2facenum, METH_VARARGS|METH_KEYWORDS, "Convert an x,y coordinate pair on the given face number to the corresponding cube face number."},
	{"radial_query_it", (PyCFunction)pyq3c_q3c_radial_query_it, METH_VARARGS|METH_KEYWORDS, ""},
	{NULL, NULL, 0, NULL}	// sentinel
};

// Ref: https://docs.python.org/3/c-api/module.html#c.PyModuleDef
static struct PyModuleDef pyq3cmodule_definition = {
	PyModuleDef_HEAD_INIT,
	"_q3c_wrapper",					// name of module
	"A Python wrapper around Q3C.",	// docstring for the module
	-1,						// size of per-interpreter state of the module,
			                // or -1 if the module keeps state in global variables.
	pyq3c_methods			// pointer to a table of module-level functions (PyMethodDef class)
//	NULL,					// m_slots
//	NULL,					// m_traverse
//	NULL,					// m_clear, a clear function to call during GC clearing of the module object
//	NULL					// m_free, function to call during deallocation of the module object
};

PyMODINIT_FUNC
PyInit__q3c_wrapper(void)
{
	return PyModule_Create(&pyq3cmodule_definition);
}


