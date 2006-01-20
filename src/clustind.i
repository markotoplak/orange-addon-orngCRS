%module orngCRS
%{
#include "interface.h"
#include "svm.h"
#include "nb.h"
#include "kikuchi.h"
%}

%typemap(python,ignore) struct LRInfo *OutValue {
    $1 = (struct LRInfo *)malloc(sizeof(struct LRInfo));
}

%typemap(python,ignore) struct CHInfo *OutValue {
    $1 = (struct CHInfo *)malloc(sizeof(struct CHInfo));
}

%typemap(python,ignore) struct CMInfo *OutValue {
    $1 = (struct CMInfo *)malloc(sizeof(struct CMInfo));
}

%typemap(python,ignore) struct CFInfo *OutValue {
    $1 = (struct CFInfo *)malloc(sizeof(struct CFInfo));
}

%typemap(python,argout) struct LRInfo *OutValue {
       int i,j;
       PyObject *o, *chisq, *devnce, *ndf, *beta, *se_beta, *fit, *covbeta, *deps; 
	   PyObject *error, *stdres, *ll, *kk;
//		printf("X+");
//		printf("lr:+1\n");
       o = PyTuple_New(10);
       chisq = PyFloat_FromDouble($1->chisq);
       devnce= PyFloat_FromDouble($1->devnce);
       ndf = PyInt_FromLong($1->ndf);
	   error = PyInt_FromLong($1->error);
       beta = PyList_New($1->k+1);
       for (i = 0; i <= $1->k; ++i) {
          ll = PyFloat_FromDouble($1->beta[i]);
		  PyList_SetItem(beta, i, ll);
       }
       se_beta = PyList_New($1->k+1);
       for (i = 0; i <= $1->k; ++i) {
          ll = PyFloat_FromDouble($1->se_beta[i]);
		  PyList_SetItem(se_beta, i, ll);
       }
       fit = PyList_New($1->nn);
       for (i = 1; i <= $1->nn; ++i) {
          ll = PyFloat_FromDouble($1->fit[i]);
		  PyList_SetItem(fit, i-1, ll);
       }
       stdres = PyList_New($1->nn);
       for (i = 1; i <= $1->nn; ++i) {
          ll = PyFloat_FromDouble($1->stdres[i]);
		  PyList_SetItem(stdres, i-1, ll);
       }
       covbeta = PyList_New($1->k+1);
       for (i = 0; i <= $1->k; ++i) {
	      kk =  PyList_New($1->k+1);
	      for(j = 0; j <= $1->k; ++j) {
			ll = PyFloat_FromDouble($1->cov_beta[i][j]);
			PyList_SetItem(kk, j, ll);
		  }
		  PyList_SetItem(covbeta, i, kk);
       }
       deps = PyList_New($1->k);
       for (i = 1; i <= $1->k; ++i) {
          ll = PyInt_FromLong($1->dependent[i]);
		  PyList_SetItem(deps, i-1, ll);
       }

       PyTuple_SetItem(o, 0, chisq);
       PyTuple_SetItem(o, 1, devnce);
       PyTuple_SetItem(o, 2, ndf);
       PyTuple_SetItem(o, 3, beta);
       PyTuple_SetItem(o, 4, se_beta);
	   PyTuple_SetItem(o, 5, fit);
       PyTuple_SetItem(o, 6, covbeta);
       PyTuple_SetItem(o, 7, stdres);
	   PyTuple_SetItem(o, 8, error);
	   PyTuple_SetItem(o, 9, deps);

	   LRInfoCleanup($1);

       if ((!$result) || ($result == Py_None)) {
	       $result = o;
       } else {
	       if (!PyList_Check($result)) {
		       PyObject *o2 = $result;
		       $result = PyList_New(0);
		       PyList_Append($result,o2);
		       Py_XDECREF(o2);
	       }
	       PyList_Append($result,o);
	       Py_XDECREF(o);
       }
//	   printf("lr:+2\n");
//		printf("-X");
}


%typemap(python,argout) struct CHInfo *OutValue {
       int i;
       PyObject *o, *ac, *n, *height, *order, *merging, *ll;

       o = PyTuple_New(5);
       n = PyInt_FromLong($1->n);
       ac = PyFloat_FromDouble($1->ac);
       merging = PyList_New(($1->n-1)*2);
       for (i = 0; i < ($1->n-1)*2; ++i) {
          ll = PyInt_FromLong($1->merging[i]);
		  PyList_SetItem(merging, i, ll);
       }
       order = PyList_New($1->n);
       for (i = 0; i < $1->n; ++i) {
          ll = PyInt_FromLong($1->order[i]);
		  PyList_SetItem(order, i, ll);
       }
       height = PyList_New($1->n-1);
       for (i = 0; i < $1->n-1; ++i) {
          ll = PyFloat_FromDouble($1->height[i+1]);
		  PyList_SetItem(height, i, ll);
       }
       PyTuple_SetItem(o, 0, n);
       PyTuple_SetItem(o, 1, merging);
       PyTuple_SetItem(o, 2, order);
       PyTuple_SetItem(o, 3, height);
       PyTuple_SetItem(o, 4, ac);

	   free($1->merging);
	   free($1->order);
	   free($1->height);
	   free($1);
       if ((!$result) || ($result == Py_None)) {
	       $result = o;
       } else {
	       if (!PyList_Check($result)) {
		       PyObject *o2 = $result;
		       $result = PyList_New(0);
		       PyList_Append($result,o2);
		       Py_XDECREF(o2);
	       }
	       PyList_Append($result,o);
	       Py_XDECREF(o);
       }
}

%typemap(python,argout) struct CMInfo *OutValue {
       int i;
       PyObject *o, *k, *cdisp, *disp, *n, *mapping, *medoids, *ll;

       o = PyTuple_New(6);
	   
       n = PyInt_FromLong($1->n);
       k = PyInt_FromLong($1->k);
       disp = PyFloat_FromDouble($1->disp);
	   
       mapping = PyList_New($1->n);
       for (i = 0; i < $1->n; ++i) {
          ll = PyInt_FromLong($1->mapping[i]);
		  PyList_SetItem(mapping, i, ll);
       }
       medoids = PyList_New($1->k);
       for (i = 0; i < $1->k; ++i) {
          ll = PyInt_FromLong($1->medoids[i]);
		  PyList_SetItem(medoids, i, ll);
       }
       cdisp = PyList_New($1->k);
       for (i = 0; i < $1->k; ++i) {
          ll = PyFloat_FromDouble($1->cdisp[i]);
		  PyList_SetItem(cdisp, i, ll);
       }
       PyTuple_SetItem(o, 0, n);
       PyTuple_SetItem(o, 1, k);
       PyTuple_SetItem(o, 2, mapping);
       PyTuple_SetItem(o, 3, medoids);
       PyTuple_SetItem(o, 4, cdisp);
       PyTuple_SetItem(o, 5, disp);
	   free($1->medoids);
	   free($1->cdisp);
	   free($1->mapping);
	   free($1);

       if ((!$result) || ($result == Py_None)) {
	       $result = o;
       } else {
	       if (!PyList_Check($result)) {
		       PyObject *o2 = $result;
		       $result = PyList_New(0);
		       PyList_Append($result,o2);
		       Py_XDECREF(o2);
	       }
	       PyList_Append($result,o);
	       Py_XDECREF(o);
       }
}


%typemap(python,argout) struct CFInfo *OutValue {
       int i,j,c;
       PyObject *o, *k, *cdisp, *disp, *n, *iterations, *value, *mapping, *membership, *ll, *kk;

       o = PyTuple_New(8);
	   
       n = PyInt_FromLong($1->n);
       k = PyInt_FromLong($1->k);
	   iterations = PyInt_FromLong($1->iterations);
	   disp = PyFloat_FromDouble($1->disp);
	   value = PyFloat_FromDouble($1->value);
	   
       mapping = PyList_New($1->n);
       for (i = 0; i < $1->n; ++i) {
          ll = PyInt_FromLong($1->mapping[i]);
		  PyList_SetItem(mapping, i, ll);
       }
       cdisp = PyList_New($1->k);
       for (i = 0; i < $1->k; ++i) {
          ll = PyFloat_FromDouble($1->cdisp[i]);
		  PyList_SetItem(cdisp, i, ll);
       }
       membership = PyList_New($1->k);
	   c = 0;
       for (i = 0; i < $1->k; ++i) {
		  kk = PyList_New($1->n);
		  for (j = 0; j < $1->n; ++j) {
			ll = PyFloat_FromDouble($1->membership[c++]);
			PyList_SetItem(kk, j, ll);
		  }
	   	  PyList_SetItem(membership, i, kk);
       }
       PyTuple_SetItem(o, 0, n);
       PyTuple_SetItem(o, 1, k);
       PyTuple_SetItem(o, 2, value);
       PyTuple_SetItem(o, 3, iterations);
       PyTuple_SetItem(o, 4, membership);
       PyTuple_SetItem(o, 5, mapping);
       PyTuple_SetItem(o, 6, cdisp);
       PyTuple_SetItem(o, 7, disp);
	   free($1->membership);
	   free($1->cdisp);
	   free($1->mapping);
	   free($1);

       if ((!$result) || ($result == Py_None)) {
	       $result = o;
       } else {
	       if (!PyList_Check($result)) {
		       PyObject *o2 = $result;
		       $result = PyList_New(0);
		       PyList_Append($result,o2);
		       Py_XDECREF(o2);
	       }
	       PyList_Append($result,o);
	       Py_XDECREF(o);
       }
}

// This tells SWIG to treat double * as a special case
%typemap(python,in) double * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int size = PyList_Size($input);
   int i = 0;
   $1 = (double *) malloc((1+size)*sizeof(double));
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyFloat_Check(o))
       $1[i] = PyFloat_AsDouble(PyList_GetItem($input,i));
     else if (PyInt_Check(o)) {
	   $1[i] = PyInt_AsLong(PyList_GetItem($input,i));
	 } else {
       PyErr_SetString(PyExc_TypeError,"list must contain floats");
       free($1);
       return NULL;
     }
   }
   $1[i] = 3.1459001;
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}

// This tells SWIG to treat int * as a special case
%typemap(python,in) int * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int size = PyList_Size($input);
   int i = 0;
   $1 = (int *) malloc((1+size)*sizeof(int));
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyInt_Check(o)) {
	   $1[i] = PyInt_AsLong(PyList_GetItem($input,i));
	 } else {
       PyErr_SetString(PyExc_TypeError,"list must contain floats");
       free($1);
       return NULL;
     }
   }
   $1[i] = -1;
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}


%typemap(python,in) struct SVMInput * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int i, j, k, size, msize;

   size = PyList_Size($input);
   msize = 0;

   $1 = (struct SVMInput *) malloc(sizeof(struct SVMInput));
   $1->data = NULL;
   $1->masking = NULL;
   $1->label = NULL;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o);
		 if ($1->data == NULL) {
			 msize = zsize;
			 $1->nn = size;
			 $1->k = msize-1;
			 $1->total = ($1->nn)*(msize-1+1); // size correction, last is class
			 $1->masking = (char **)malloc(sizeof(char *)*(size+1));
			 for (j = 0; j <= size; ++j) {
				$1->masking[j] = (char *)malloc(sizeof(char)*($1->k+1));
				for (k = 0; k <= $1->k; ++k) {
					$1->masking[j][k] = 0;
				}
			 }
			 $1->data = (double **)malloc(sizeof(double *)*(size+1));
			 for (j = 0; j <= size; ++j) {
				$1->data[j] = (double *)malloc(sizeof(double)*($1->k+1));
			 }
			 $1->label = (double *)malloc(sizeof(double)*($1->nn+1));
		 }
		 if (zsize == msize) {
			 // fetch the class (must be binary)
			 PyObject *p = PyList_GetItem(o,zsize-1);
			 if (p == NULL) {
				 PyErr_SetString(PyExc_TypeError,"cannot fetch example");
				 SVMcleanup($1);
				 return NULL;
			 }

			 if (PyInt_Check(p) || PyLong_Check(p)) {
				$1->label[i] = (double)PyInt_AsLong(p);
			 } else if(PyFloat_Check(p)) {
				$1->label[i] = PyFloat_AsDouble(p);
			 } else {
				 PyErr_SetString(PyExc_TypeError,"class not float or int");
				 SVMcleanup($1);
				 return NULL;
			 }
				
			 // fetch the attribute values
			 for (j = 0; j < zsize-1; ++j) {
			     PyObject *p = PyList_GetItem(o,j);

	 			 if(PyTuple_Check(p)) {
					 PyObject *ps = PyTuple_GetItem(p,1); // masking
					 if (PyInt_Check(ps) || PyLong_Check(ps)) {
						$1->masking[i][j] = (char)PyInt_AsLong(ps);
					 } else if(PyFloat_Check(ps)) {
						$1->masking[i][j] = (char)PyFloat_AsDouble(ps);
					 } else {
						 PyErr_SetString(PyExc_TypeError,"masking not 1 or 0");
						 SVMcleanup($1);
						 return NULL;
					 }
					 p = PyTuple_GetItem(p,0); // foist value
				 }
				 if (PyFloat_Check(p)) {
					// correct +1 because that array starts at 1
					$1->data[i][j] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[i][j] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"examples must contain doubles or ints. if masking, use tuples");
					 SVMcleanup($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"examples must be of equal size");
			 SVMcleanup($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"example table must contain examples as lists");
       SVMcleanup($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
 //("svmi:done\n");
}

%typemap(python,in) struct SVMExample * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int j, msize;

   msize = PyList_Size($input);

   $1 = (struct SVMExample *)malloc(sizeof(struct SVMExample));
   $1->k = msize;
   $1->data    = (double *)malloc(sizeof(double)*msize);
   $1->masking = (char *)malloc(sizeof(char)*msize);

	// fetch the attribute values
	for (j = 0; j < msize; ++j) {
	 PyObject *p = PyList_GetItem($input,j);

     $1->masking[j] = 0;
	 if(PyTuple_Check(p)) {
		$1->masking[j] = 1;
		 p = PyTuple_GetItem(p,0); // foist value
	 }
	 if (PyFloat_Check(p)) {
		$1->data[j] = PyFloat_AsDouble(p);
	 } else if (PyInt_Check(p)) {
		$1->data[j] = (double)PyInt_AsLong(p);
	 } else {
		 PyErr_SetString(PyExc_TypeError,"examples must contain doubles or ints. if masking, use tuples");
		 SVMexcleanup($1);
		 return NULL;
	 }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}


%typemap(python,in) struct SVMSparseInput * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int i, j, k, nex, nvals;
   
   nex = PyList_Size($input);
   nvals = 0;

   $1 = (struct SVMSparseInput *) malloc(sizeof(struct SVMSparseInput));
   $1->elements = 0;
   $1->index = NULL;
   $1->value = NULL;
   $1->label = NULL;
   for (i = 0; i < nex; i++) {
     PyObject *ex = PyList_GetItem($input,i);
     // example
     if (PyList_Check(ex)) {
		 nvals = PyList_Size(ex)-1;
		 if ($1->index == NULL) {
			 $1->nn = nex;
			 $1->lengths = (int *)malloc($1->nn*sizeof(int));
			 $1->value = (double **)malloc($1->nn*sizeof(double *));
			 $1->index = (int **)malloc($1->nn*sizeof(int *));
			 $1->label = (double *)malloc($1->nn*sizeof(double));
			 for(j = 0; j < nex; ++j)
				$1->lengths[j] = -1;
		 }
		// fetch the class (must be binary)
		PyObject *p = PyList_GetItem(ex,0);
		if (p == NULL) {
			PyErr_SetString(PyExc_TypeError,"cannot fetch label");
			SVMscleanup($1);
 			return NULL;
 		}
		if (PyInt_Check(p) || PyLong_Check(p)) {
			$1->label[i] = (double)PyInt_AsLong(p);
		} else if(PyFloat_Check(p)) {
			$1->label[i] = PyFloat_AsDouble(p);
		} else {
			PyErr_SetString(PyExc_TypeError,"class not float or int");
			SVMscleanup($1);
			return NULL;
		}
		$1->lengths[i] = nvals-1;
		$1->elements += nvals-1;
		if(nvals > 0) {
			$1->value[i] = (double *)malloc((nvals-1)*sizeof(double));
			$1->index[i] = (int *)malloc((nvals-1)*sizeof(int));
		}
		// fetch the attribute values
		for (j = 1; j < nvals; ++j) {
			PyObject *p = PyList_GetItem(ex,j);
	 		if(PyTuple_Check(p)) {
	 			PyObject *idx = PyTuple_GetItem(p,0); 
				PyObject *val = PyTuple_GetItem(p,1); 
				if ((PyInt_Check(idx) || PyLong_Check(idx)) && (PyInt_AsLong(idx) >= 1) ) {
					$1->index[i][j-1] = PyInt_AsLong(idx);
				} else {
					PyErr_SetString(PyExc_TypeError,"index not a positive integer");
					SVMscleanup($1);
					return NULL;
				}
				if (PyInt_Check(val) || PyLong_Check(val)) {
					$1->value[i][j-1] = (double)PyInt_AsLong(val);
				} else if(PyFloat_Check(val)) {
					$1->value[i][j-1] = PyFloat_AsDouble(val);
				}  else {
					PyErr_SetString(PyExc_TypeError,"value not an integer or a float");
					SVMscleanup($1);
					return NULL;
				}
			} else {
				PyErr_SetString(PyExc_TypeError,"attribute values within an example must be tuples (index,value)");
				SVMscleanup($1);
				return NULL;
			}
		 } 
     } else {
       PyErr_SetString(PyExc_TypeError,"example table must contain examples as lists of tuples");
       SVMscleanup($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}

%typemap(python,in) struct SVMSparseExample * {
int nvals, i,j; 
 /* Check if is a list */
     if (PyList_Check($input)) {
		 nvals = PyList_Size($input);
		 $1 = (SVMSparseExample *)malloc(sizeof(SVMSparseExample));
		 $1->nn = nvals;
 		 $1->value = (double *)malloc((nvals)*sizeof(double));
		 $1->index = (int *)malloc((nvals)*sizeof(int));
		 // fetch the attribute values
		 for (j = 0; j < nvals; ++j) {
			PyObject *p = PyList_GetItem($input,j);
	 		if(PyTuple_Check(p)) {
	 			PyObject *idx = PyTuple_GetItem(p,0); 
				PyObject *val = PyTuple_GetItem(p,1); 
				if (PyInt_Check(idx) || PyLong_Check(idx)) {
					$1->index[j] = PyInt_AsLong(idx);
				} else {
					PyErr_SetString(PyExc_TypeError,"index not 1 or 0");
					SVMsecleanup($1);
					return NULL;
				}
				if (PyInt_Check(val) || PyLong_Check(val)) {
					$1->value[j] = (double)PyInt_AsLong(val);
				} else if(PyFloat_Check(val)) {
					$1->value[j] = PyFloat_AsDouble(val);
				}  else {
					PyErr_SetString(PyExc_TypeError,"value not an integer or a float");
					SVMsecleanup($1);
					return NULL;
				}
			} else { 
				PyErr_SetString(PyExc_TypeError,"attribute values within an example must be tuples (index,value)");
				SVMsecleanup($1);
				return NULL;
			}
		 } 
     } else {
       PyErr_SetString(PyExc_TypeError,"example must be a list of tuples");
       SVMsecleanup($1);
       return NULL;
     }
}


%typemap(python,in) struct LRInput * {
 /* Check if is a list */
// printf("A-");
 if (PyList_Check($input)) {
   int i, j, size, msize;

   size = PyList_Size($input);
   msize = 0;
	//printf("lr:-1\n");
   $1 = (struct LRInput *) malloc(sizeof(struct LRInput));
   $1->data = NULL;
   $1->success = NULL;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o);
		 if ($1->data == NULL) {
			 msize = zsize;
			 $1->nn = size;
			 $1->k = msize-1; // size correction, last is class
			 $1->data = (double **)malloc(sizeof(double *)*(size+1));
			 for (j = 0; j <= size; ++j) {
				$1->data[j] = (double *)malloc(sizeof(double)*($1->k+1));
			 }
			 $1->success = (double *)malloc(sizeof(double)*($1->nn+1));
			 $1->trials = (double *)malloc(sizeof(double)*($1->nn+1));
		 }
		 if (zsize == msize) {
			 // fetch the class (must be either binary or a tuple)
			 PyObject *p = PyList_GetItem(o,zsize-1);
			 if (p == NULL) {
				 PyErr_SetString(PyExc_TypeError,"cannot fetch example");
				 LRcleanup($1);
				 return NULL;
			 }
				
			 // the way of storing the class value is either simple scalar 
			 // probability, or tuple (probability, weight)

			 // set the default value for # trials
			 $1->trials[i+1] = 1.0;
	 		 if(PyTuple_Check(p)) {
				 PyObject *ps = PyTuple_GetItem(p,1); // get # trials
				 if (PyInt_Check(ps) || PyLong_Check(ps)) {
					$1->trials[i+1] = (double)PyInt_AsLong(ps);
				 } else if(PyFloat_Check(ps)) {
					$1->trials[i+1] = PyFloat_AsDouble(ps);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"number of trials is not float or int");
					 LRcleanup($1);
					 return NULL;
				 }
				 p = PyTuple_GetItem(p,0); // foist a value for successes
			 } 
				
			 if (PyInt_Check(p) || PyLong_Check(p)) {
				$1->success[i+1] = $1->trials[i+1]*(double)PyInt_AsLong(p);
			 } else if(PyFloat_Check(p)) {
				$1->success[i+1] = $1->trials[i+1]*PyFloat_AsDouble(p);
			 } else {
				 PyErr_SetString(PyExc_TypeError,"experiment success is not float or int");
				 LRcleanup($1);
				 return NULL;
			 }
			 				
			 // fetch the attribute values
			 for (j = 0; j < zsize-1; ++j) {
			     PyObject *p = PyList_GetItem(o,j);
				 if (PyFloat_Check(p)) {
					// correct +1 because that array starts at 1
					$1->data[i+1][j+1] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[i+1][j+1] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"examples must contain doubles or ints");
					 LRcleanup($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"examples must be of equal size");
			 LRcleanup($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"example table must contain examples as lists");
       LRcleanup($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
//  printf("-A");
 //printf("lr:-2\n");
}

%typemap(python,in) struct CInput * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int i, j, size, msize;

   size = PyList_Size($input);
   msize = 0;

   $1 = (struct CInput *) malloc(sizeof(struct CInput));
   $1->data = NULL;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o);
		 if ($1->data == NULL) {
			 msize = zsize;
			 $1->nn = size;
			 $1->jpp = msize;
			 $1->data = (double *)malloc(sizeof(double)*zsize*size);
		 }
		 if (zsize == msize) {
			 for (j = 0; j < zsize; ++j) {
			     PyObject *p = PyList_GetItem(o,j);
				 if (PyFloat_Check(p)) {
					$1->data[i+j*size] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[i+j*size] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"sublists must contain doubles");
					 if ($1->data != NULL)
						 free($1->data);
					 free($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"sublists must be of equal size");
			 if ($1->data != NULL)
				 free($1->data);
			 free($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"list must contain lists");
	   if ($1->data != NULL)
		   free($1->data);
       free($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}


%typemap(python,in) struct DInput * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int i, j, size, offset;

   size = PyList_Size($input);     /* the number of elements - 1 */

   $1 = (struct DInput *) malloc(sizeof(struct DInput));
   $1->nn = size+1;
   $1->data = (double *)malloc((1 + (size * (size + 1))/2)*sizeof(double));

   $1->data[0] = 0.0; // offset for pam and agnes
   //offset = (size * (size + 1))/2;
   offset = 1;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o); /* the line of a matrix */
		 if (zsize == i+1) {
			 for (j = 0; j < zsize; ++j) {
			     PyObject *p = PyList_GetItem(o,j);
				 if (PyFloat_Check(p)) {
					$1->data[offset++] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[offset++] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"sublists must contain doubles");
					 if ($1->data != NULL)
						 free($1->data);
					 free($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"list must contain the bottom-triangular dissimilarity matrix");
			 if ($1->data != NULL)
				 free($1->data);
			 free($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"list must contain lists");
	   if ($1->data != NULL)
		   free($1->data);
       free($1);
       return NULL;
     }
   }
   //assert(offset == 0);
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}


%typemap(python,arginit) double *{
  $1 = NULL;
}
%typemap(python,arginit) int *{
  $1 = NULL;
}
%typemap(python,arginit) struct LRInput *{
  $1 = NULL;
}
%typemap(python,arginit) struct SVMInput *{
  $1 = NULL;
}
%typemap(python,arginit) struct SVMExample * {
  $1 = NULL;
}
%typemap(python,arginit) struct SVMSparseInput *{
  $1 = NULL;
}
%typemap(python,arginit) struct SVMSparseExample * {
  $1 = NULL;
}
%typemap(python,arginit) struct CInput * {
  $1 = NULL;
}
%typemap(python,arginit) struct DInput * {
  $1 = NULL;
}

// This cleans up the double * array we malloc'd before the function call
%typemap(python,freearg) double * {
  if ($1 != NULL)
	free((double *) $1);
}

// This cleans up the double * array we malloc'd before the function call
%typemap(python,freearg) int * {
  if ($1 != NULL)
	  free((int *) $1);
}

%typemap(python,freearg) struct LRInput * {
  if($1 != NULL) {
	  int i;
	  for (i = 0; i <= $1->nn; ++i)
		free((double*)$1->data[i]);
	  free((double**) $1->data);
	  free((double*) $1->success);
	  free((double*) $1->trials);
	  free((struct LRInput*) $1);
  }
}

%typemap(python,freearg) struct SVMInput * {
  if($1 != NULL) {
	  int i;
	  //printf("svmc:1\n");
	  for (i = 0; i <= $1->nn; ++i) {
		free((double*)$1->data[i]);
		free((char*)$1->masking[i]);
	  }
	  free((double**) $1->data);
	  free((char**) $1->masking);
	  free((double*) $1->label);
	  free((struct SVMInput*) $1);
	  //("svmc:2\n");
  }
}

%typemap(python,freearg) struct SVMExample * {
  if ($1 != NULL) {
	  free((double*) $1->data);
	  free((char*) $1->masking);
	  free((struct SVMExample*) $1);
  }
}

%typemap(python,freearg) struct SVMSparseInput * {
  if($1 != NULL) {
	  int i;
	  for (i = 0; i < $1->nn; ++i) {
	    if($1->lengths[i] > 0) {
			free((double*)$1->value[i]);
			free((int*)$1->index[i]);
		}
	  }
	  free((double**) $1->value);
	  free((int**) $1->index);
	  free((double*) $1->label);
	  free((int*) $1->lengths);
	  free((struct SVMSparseInput*) $1);
  }
}

%typemap(python,freearg) struct SVMSparseExample * {
  if ($1 != NULL) {
	  free((double*) $1->value);
	  free((int *) $1->index);
	  free((struct SVMSparseExample*) $1);
  }
}

%typemap(python,freearg) struct CInput * {
  if ($1 != NULL) {
	  free((double*) $1->data);
	  free((struct CInput*) $1);
  }
}

%typemap(python,freearg) struct DInput * {
  if ($1 != NULL) {
	  free((double*) $1->data);
	  free((struct DInput*) $1);
  }
}


%typemap(python, in) struct svm_model * {
    PyObject *o, *t, *p, *zz, *uu;
	int m,l,i,j, elements, pairs;
	struct svm_node *x_space;
	struct svm_model *model;
	struct svm_parameter *param;


	o = $input;
	model = (struct svm_model *)malloc(sizeof(struct svm_model));
	param = &(model->param);
	model->label = NULL;
	model->nSV = NULL;
	model->probA = NULL;
	model->probB = NULL;

	if (!PyDict_Check(o)) {
       PyErr_SetString(PyExc_TypeError,"must be a dictionary");
       free(model);
       return NULL;
    }
	t = PyDict_GetItemString(o, "svm_type");
	if (!PyInt_Check(t)) {
       PyErr_SetString(PyExc_TypeError,"svm_type missing");
       free(model);
       return NULL;
    }
	param->svm_type = PyInt_AsLong(t);
	t = PyDict_GetItemString(o, "kernel_type");
	if (!PyInt_Check(t)) {
       PyErr_SetString(PyExc_TypeError,"kernel_type missing");
       free(model);
       return NULL;
    }
	param->kernel_type = PyInt_AsLong(t);
    	
	if(param->kernel_type == POLY) {
		t = PyDict_GetItemString(o, "degree");
		if (!PyFloat_Check(t)) {
		   PyErr_SetString(PyExc_TypeError,"degree missing");
		   free(model);
		   return NULL;
		}
		param->degree = PyFloat_AsDouble(t);
	}
	
	if(param->kernel_type == POLY || param->kernel_type == RBF || param->kernel_type == SIGMOID) {
		t = PyDict_GetItemString(o, "gamma");
		if (!PyFloat_Check(t)) {
		   PyErr_SetString(PyExc_TypeError,"gamma missing");
		   free(model);
		   return NULL;
		}
		param->gamma = PyFloat_AsDouble(t);
	}

	if(param->kernel_type == POLY || param->kernel_type == SIGMOID) {
		t = PyDict_GetItemString(o, "coef0");
		if (!PyFloat_Check(t)) {
		   PyErr_SetString(PyExc_TypeError,"coef0 missing");
		   free(model);
		   return NULL;
		}
		param->coef0 = PyFloat_AsDouble(t);
	}

	t = PyDict_GetItemString(o, "nr_class");
	if (!PyInt_Check(t)) {
	   PyErr_SetString(PyExc_TypeError,"nr_class missing");
	   free(model);
	   return NULL;
	}
	m = model->nr_class = PyInt_AsLong(t);
	t = PyDict_GetItemString(o, "total_sv");
	if (!PyInt_Check(t)) {
	   PyErr_SetString(PyExc_TypeError,"total_sv missing");
	   free(model);
	   return NULL;
	}
	l = model->l = PyInt_AsLong(t);


	t = PyDict_GetItemString(o, "rho");
	if (!PyList_Check(t)) {
	   PyErr_SetString(PyExc_TypeError,"rho missing");
	   free(model);
	   return NULL;
	}
	pairs = m*(m-1)/2;
	model->rho = (double *)malloc(sizeof(double)*pairs);
	for(i=0;i<pairs;i++) {
		p = PyList_GetItem(t, i);
		model->rho[i] = PyFloat_AsDouble(p);
	}

	t = PyDict_GetItemString(o, "ProbA");
	if (t != NULL) {
		model->probA = (double *)malloc(sizeof(double)*pairs);
		for(i=0;i<pairs;i++) {
			p = PyList_GetItem(t, i);
			model->probA[i] = PyFloat_AsDouble(p);
		}
	}

	t = PyDict_GetItemString(o, "ProbB");
	if (t != NULL) {
		model->probB = (double *)malloc(sizeof(double)*pairs);
		for(i=0;i<pairs;i++) {
			p = PyList_GetItem(t, i);
			model->probB[i] = PyFloat_AsDouble(p);
		}
	}

	t = PyDict_GetItemString(o, "label");
	if (t != NULL) {
		model->label = (int *)malloc(sizeof(int)*m);
		for(i=0;i<m;i++) {
			p = PyList_GetItem(t, i);
			model->label[i] = PyInt_AsLong(p);
		}
	}

	t = PyDict_GetItemString(o, "nr_sv");
	if (t != NULL) {
		model->nSV = (int *)malloc(sizeof(int)*m);
		for(i=0;i<m;i++) {
			p = PyList_GetItem(t, i);
			model->nSV[i] = PyInt_AsLong(p);
		}
	}

	t = PyDict_GetItemString(o, "elements");
	elements = PyInt_AsLong(t);
	
	model->sv_coef = (double **)malloc(sizeof(double *)*m);
	model->SVidx = (int *)malloc(sizeof(int)*l); // allocate even if not use
	model->SV = (struct svm_node **)malloc(sizeof(struct svm_node **)*l);
	for(i=0;i<m;i++)
		model->sv_coef[i] = (double *)malloc(sizeof(double)*l);

	x_space = (struct svm_node *)malloc(sizeof(struct svm_node)*elements);

	p = PyDict_GetItemString(o, "SV");
	if (!PyList_Check(p) || PyList_Size(p)!=l) {
	   PyErr_SetString(PyExc_TypeError,"SV list missing");
	   free(model);
	   return NULL;
	}
	j = 0;
	for(i=0;i<l;i++)
	{
		int jj, kk;

		zz = PyList_GetItem(p,i);
		if (!PyList_Check(zz)) {
		   PyErr_SetString(PyExc_TypeError,"wrong SV vector (leak)"); return NULL;
		}

		t = PyList_GetItem(zz, 0); // sv_coef is first
		if (!PyList_Check(t) || PyList_Size(t)!= m-1 ) {
		   PyErr_SetString(PyExc_TypeError,"SV coef wrong (leak)"); return NULL;
		}

		for(jj=0;jj<m-1;jj++) {
			uu = PyList_GetItem(t, jj);
			if (!PyFloat_Check(uu)) {
				PyErr_SetString(PyExc_TypeError,"SV coef entry wrong (leak)"); return NULL;
			}
			model->sv_coef[jj][i] = PyFloat_AsDouble(uu);
		}

		model->SV[i] = &(x_space[j]);

		kk = PyList_Size(zz);

		for (jj = 1; jj < kk; ++jj) {
			t =PyList_GetItem(zz, jj);
			if (!PyTuple_Check(t)) {
				PyErr_SetString(PyExc_TypeError,"SV entry wrong (leak)"); return NULL;
			}
			
			uu = PyTuple_GetItem(t,0);
			x_space[j].index = PyInt_AsLong(uu);
			uu = PyTuple_GetItem(t,1);
			x_space[j].value = PyFloat_AsDouble(uu);

			++j;
		}
		x_space[j++].index = -1;
	}
	assert(j == elements);
	model->free_sv = 1;	// XXX
	$1 = model;
}

%typemap(python, out) struct wsvm_model * {
	struct svm_parameter *param = &($1->m->param);
	double **sv_coef;
	struct svm_node **SV;
	struct svm_model *model = $1->m;
    int i, j, nr_class, l, elements;
    PyObject *o, *t, *p, *ip;

	o = PyDict_New();
    t = PyInt_FromLong(param->svm_type);
	PyDict_SetItemString(o, "svm_type", t); Py_XDECREF(t);
    t = PyInt_FromLong(param->kernel_type);
	PyDict_SetItemString(o, "kernel_type", t); Py_XDECREF(t);
    	
	if(param->kernel_type == POLY) {
		t =	PyFloat_FromDouble(param->degree);
		PyDict_SetItemString(o, "degree", t); Py_XDECREF(t);
	}

	if(param->kernel_type == POLY || param->kernel_type == RBF || param->kernel_type == SIGMOID) {
		t = PyFloat_FromDouble(param->gamma);
		PyDict_SetItemString(o, "gamma", t); Py_XDECREF(t);
	}

	if(param->kernel_type == POLY || param->kernel_type == SIGMOID) {
		t = PyFloat_FromDouble(param->coef0);
		PyDict_SetItemString(o, "coef0", t); Py_XDECREF(t);
	}

	nr_class = model->nr_class;
	l = model->l;

    t = PyInt_FromLong(nr_class);
	PyDict_SetItemString(o, "nr_class", t); Py_XDECREF(t);
    t = PyInt_FromLong(l);
	PyDict_SetItemString(o, "total_sv", t); Py_XDECREF(t);

	{
	   	t = PyList_New(nr_class*(nr_class-1)/2);
		for(i=0;i<nr_class*(nr_class-1)/2;i++) {
			PyList_SetItem(t, i, PyFloat_FromDouble(model->rho[i]));
		}
		PyDict_SetItemString(o, "rho", t); Py_XDECREF(t);
	}
	
	if(model->probA) {
	   	t = PyList_New(nr_class*(nr_class-1)/2);
		for(i=0;i<nr_class*(nr_class-1)/2;i++) {
			PyList_SetItem(t, i, PyFloat_FromDouble(model->probA[i]));
		}
		PyDict_SetItemString(o, "ProbA", t); Py_XDECREF(t);
	}
	
	if(model->probB) {
	   	t = PyList_New(nr_class*(nr_class-1)/2);
		for(i=0;i<nr_class*(nr_class-1)/2;i++) {
			PyList_SetItem(t, i, PyFloat_FromDouble(model->probB[i]));
		}
		PyDict_SetItemString(o, "ProbB", t); Py_XDECREF(t);
	}
	
	if(model->label) {
	   	t = PyList_New(nr_class);
		for(i=0;i<nr_class;i++) {
			PyList_SetItem(t, i, PyInt_FromLong(model->label[i]));
		}
		PyDict_SetItemString(o, "label", t); Py_XDECREF(t);
	}
	
	
	if(model->nSV) {
	   	t = PyList_New(nr_class);
		for(i=0;i<nr_class;i++) {
			PyList_SetItem(t, i, PyInt_FromLong(model->nSV[i]));
		}
		PyDict_SetItemString(o, "nr_sv", t); Py_XDECREF(t);
	}
	
	
	sv_coef = model->sv_coef;
	SV = model->SV;
	
	p = PyList_New(l);
	ip = PyList_New(l);
	elements = l; // terminals
	for(i=0;i<l;i++)
	{
		const struct svm_node *pp;
		PyObject *zz;

	   	t = PyList_New(nr_class-1);
		for(j=0;j<nr_class-1;j++) {
			PyList_SetItem(t, j, PyFloat_FromDouble(model->sv_coef[j][i]));
		}

		// count attributes
		pp = SV[i];
		PyList_SetItem(ip, i, PyInt_FromLong(model->SVidx[i])); // save the SV index
		j = 1;
		while(pp->index != -1)
		{
			pp++; ++j;
		}

		// allocate
		zz = PyList_New(j);
		PyList_SetItem(zz, 0, t); // sv_coef is first

		pp = SV[i];
		j = 0;
		while(pp->index != -1)
		{
			t = PyTuple_New(2);
			PyTuple_SetItem(t,0,PyInt_FromLong(pp->index));
			PyTuple_SetItem(t,1,PyFloat_FromDouble(pp->value));
			pp++; ++j; ++elements;
			PyList_SetItem(zz,j,t);
		}

		PyList_SetItem(p,i,zz);
	}
	PyDict_SetItemString(o, "SV", p); Py_XDECREF(p);
	PyDict_SetItemString(o, "SVi", ip); Py_XDECREF(ip);

//  printf("svmo:2\n");
    t = PyInt_FromLong(elements);
	PyDict_SetItemString(o, "elements", t); Py_XDECREF(t);
	svm_destroy_model($1->m);
	free($1->prob->x);	
	free($1->prob->y);	
	free($1->prob);	
	free($1->x_space);	
	free($1);
    $result = o;
} 



%typemap(python,ignore) struct SVMOut *OutValue {
    $1 = (struct SVMOut *)malloc(sizeof(struct SVMOut));
}
%typemap(python,argout) struct SVMOut *OutValue {
   int i;
   PyObject *o, *q;

   o = PyList_New($1->nn);

   for(i = 0; i < $1->nn; ++i) {
	 q = PyFloat_FromDouble($1->v[i]);
	 PyList_SetItem(o, i, q);
   }
   free($1->v);
   free($1);
   if ((!$result) || ($result == Py_None)) {
	   $result = o;
   } else {
	   if (!PyList_Check($result)) {
		   PyObject *o2 = $result;
		   $result = PyList_New(0);
		   PyList_Append($result,o2);
		   Py_XDECREF(o2);
	   }
	   PyList_Append($result,o);
	   Py_XDECREF(o);
   }
}


void svm_destroy_model(psvm_model model);
psvm_model SVMClassifier(struct svm_model *InValue);
struct wsvm_model *SVMLearnS(struct SVMSparseInput *input, int svm_type, int kernel_type, double degree,
		 double gamma, double coef0, double nu, double cache_size, double C, 
		 double eps, double p, int shrinking, int probability, int nr_weight, double *weight, 
		 int *weight_label);
struct wsvm_model *SVMLearn(struct SVMInput *input, int svm_type, int kernel_type, double degree,
		 double gamma, double coef0, double nu, double cache_size, double C, 
		 double eps, double p, int shrinking, int probability, int nr_weight, double *weight, 
		 int *weight_label);
double SVMClassify(psvm_model model, struct SVMExample *input);
void SVMClassifyP(psvm_model model, struct SVMExample *input, struct SVMOut *OutValue );
void SVMClassifyM(psvm_model model, struct SVMExample *input, struct SVMOut *OutValue );
double SVMClassifyS(psvm_model model, struct SVMSparseExample *input);
void SVMClassifyPS(psvm_model model, struct SVMSparseExample *input, struct SVMOut *OutValue );
void SVMClassifyMS(psvm_model model, struct SVMSparseExample *input, struct SVMOut *OutValue );


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


%typemap(python,ignore) struct NBResult *OutValue {
    $1 = (struct NBResult *)malloc(sizeof(struct NBResult));
}

%typemap(python,in) struct NBInput * {
 int i, j, size, msize;
 PyObject *arr,*cvl,*cards,*alpha;
 /* Check if is a list */

 arr = PyTuple_GetItem($input,0);
 if (PyList_Check(arr)) {

   size = PyList_Size(arr);
   msize = 0;
	//printf("lr:-1\n");
   $1 = (struct NBInput *) malloc(sizeof(struct NBInput));
   $1->data = NULL;
   $1->l = NULL;
   $1->card = NULL;
   $1->cvi = NULL;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem(arr,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o);
		 if ($1->data == NULL) {
			 msize = zsize;
			 $1->nn = size;
			 $1->k = msize;
			 $1->na = $1->k-1; // number of attributes
			 $1->data = (int **)malloc(sizeof(int *)*(size));
			 for (j = 0; j < size; ++j) {
				$1->data[j] = (int *)malloc(sizeof(int)*($1->k));
			 }
			 $1->l = (int *)malloc(sizeof(int)*($1->nn));
			 $1->cvi = (int *)malloc(sizeof(int)*($1->nn));
			 $1->card = (int *)malloc(sizeof(int)*($1->k));
		 }
		 if (zsize == msize) {		 				
			 // fetch the attribute values
			 for (j = 0; j < $1->k; ++j) {
			     PyObject *p = PyList_GetItem(o,j);
				 if (PyFloat_Check(p)) {
					// correct +1 because that array starts at 1
					$1->data[i][j] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[i][j] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"examples must contain doubles or ints");
					 NBcleanup($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"examples must be of equal size");
			 NBcleanup($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"example table must contain examples as lists");
       NBcleanup($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"tuple element #0 not a list");
   return NULL;
 }


 // fetch the cardinalities
 cards = PyTuple_GetItem($input,1);
 if (!PyList_Check(cards) && PyList_Size(cards) != $1->k) {
   PyErr_SetString(PyExc_TypeError,"tuple element #1 not an attribute cardinality list");
   return NULL;
 }

 for (j = 0; j < $1->k; ++j) {
	 PyObject *p = PyList_GetItem(cards,j);
	 if (PyFloat_Check(p)) {
		// correct +1 because that array starts at 1
		$1->card[j] = PyFloat_AsDouble(p);
	 } else if (PyInt_Check(p)) {
		$1->card[j] = PyInt_AsLong(p);
	 } else {
		 PyErr_SetString(PyExc_TypeError,"card list must contain doubles or ints");
		 NBcleanup($1);
		 return NULL;
	 }
 }

 // fetch the cvl 
 cvl = PyTuple_GetItem($input,2);
 if (!PyList_Check(cvl) && PyList_Size(cvl) != $1->nn) {
   PyErr_SetString(PyExc_TypeError,"tuple element #2 not a CVI list");
   return NULL;
 }

 $1->maxcvi = 0;
 for (j = 0; j < $1->nn; ++j) {
	 PyObject *p = PyList_GetItem(cvl,j);
	 if (PyFloat_Check(p)) {
		// correct +1 because that array starts at 1
		$1->cvi[j] = PyFloat_AsDouble(p);
	 } else if (PyInt_Check(p)) {
		$1->cvi[j] = PyInt_AsLong(p);
	 } else {
		 PyErr_SetString(PyExc_TypeError,"cvl must contain doubles or ints");
		 NBcleanup($1);
		 return NULL;
	 }
	 if($1->maxcvi < $1->cvi[j])
		$1->maxcvi = $1->cvi[j];
 }

 alpha = PyTuple_GetItem($input,3);
 if (PyFloat_Check(alpha)) {
	$1->a = PyFloat_AsDouble(alpha);
 } else if (PyInt_Check(alpha)) {
	$1->a = PyInt_AsLong(alpha);
 } else {
	$1->a = 1.0; // default smoothing
 }
}

%typemap(python,in) struct NBModel * {
	PyObject *p, *sing, *grup;

	$1 = (struct NBModel *) malloc(sizeof(struct NBInput));
	$1->nosingles = 0;
	$1->singles = NULL;
	$1->nogroups = 0;
	$1->groups = NULL;

	sing = PyTuple_GetItem($input,0);
	if (PyList_Check(sing)) {
		int i;

		$1->nosingles = PyList_Size(sing);
		if($1->nosingles > 0) {
			$1->singles = (int *)malloc(sizeof(int)*$1->nosingles);
			for (i = 0; i < $1->nosingles; i++) {
				$1->singles[i] = PyInt_AsLong(PyList_GetItem(sing,i));
			}
		}
	} else {
		PyErr_SetString(PyExc_TypeError,"tuple element #0 not a list of singles");
		return NULL;
	}

	grup = PyTuple_GetItem($input,1);
	if (PyList_Check(grup)) {
		int i, j, msize;

		$1->nogroups = PyList_Size(grup);
		if($1->nogroups > 0) {
			$1->groups = (struct NBGroup *)malloc(sizeof(struct NBGroup)*$1->nogroups);
			for (i = 0; i < $1->nogroups; i++) {
				PyObject *o = PyList_GetItem(grup,i);
				msize = PyList_Size(o);
				$1->groups[i].n = msize;
				$1->groups[i].l = (int *)malloc(sizeof(int)*msize);
				for (j = 0; j < msize; ++j)
					$1->groups[i].l[j] = PyInt_AsLong(PyList_GetItem(o,j));
			}
		}
	} else {
		PyErr_SetString(PyExc_TypeError,"tuple element #1 not a list of groups");
		return NULL;
	}
}

%typemap(python,argout) struct NBResult *OutValue {
   int i,j,c;
   PyObject *o, *err, *q;

   o = PyTuple_New(6);

   q = PyFloat_FromDouble($1->kl_q);
   err = PyFloat_FromDouble($1->kl_err);

   PyTuple_SetItem(o, 0, q);
   PyTuple_SetItem(o, 1, err);

   q = PyFloat_FromDouble($1->b_q);
   err = PyFloat_FromDouble($1->b_err);

   PyTuple_SetItem(o, 2, q);
   PyTuple_SetItem(o, 3, err);

   q = PyFloat_FromDouble($1->er_q);
   err = PyFloat_FromDouble($1->er_err);

   PyTuple_SetItem(o, 4, q);
   PyTuple_SetItem(o, 5, err);
   free($1);

   if ((!$result) || ($result == Py_None)) {
	   $result = o;
   } else {
	   if (!PyList_Check($result)) {
		   PyObject *o2 = $result;
		   $result = PyList_New(0);
		   PyList_Append($result,o2);
		   Py_XDECREF(o2);
	   }
	   PyList_Append($result,o);
	   Py_XDECREF(o);
   }
}


%typemap(python,ignore) struct NBList *OutValue {
    $1 = (struct NBList *)malloc(sizeof(struct NBList));
}

%typemap(python,argout) struct NBList *OutValue {
   int i;
   PyObject *o, *q;

   o = PyList_New($1->n);

   for(i = 0; i < $1->n; ++i) {
	 q = PyFloat_FromDouble($1->l[i]);
	 PyList_SetItem(o, i, q);
   }
   free($1->l);
   free($1);

   if ((!$result) || ($result == Py_None)) {
	   $result = o;
   } else {
	   if (!PyList_Check($result)) {
		   PyObject *o2 = $result;
		   $result = PyList_New(0);
		   PyList_Append($result,o2);
		   Py_XDECREF(o2);
	   }
	   PyList_Append($result,o);
	   Py_XDECREF(o);
   }
}


%typemap(python,arginit) struct NBModel * {
	$1 = NULL;
}

%typemap(python,freearg) struct NBModel * {
	int i;
	if($1->singles != NULL)
	  free($1->singles);
	if($1->groups != NULL) {
	  for (i = 0; i < $1->nogroups; ++i) {
		free($1->groups[i].l);
	  }
	  free($1->groups);
	}
	free($1);
}



void MCluster(struct CInput *in, long k, int metric, struct CMInfo *OutValue);
void HCluster(struct CInput *in, int metric, int method, struct CHInfo *OutValue);
void DMCluster(struct DInput *in, long k, struct CMInfo *OutValue);
void DHCluster(struct DInput *in, int method, struct CHInfo *OutValue);
void FCluster(struct CInput *in, long k, int metric, struct CFInfo *OutValue);
void DFCluster(struct DInput *in, long k, struct CFInfo *OutValue);

void LogReg(struct LRInput *input, double regularization, struct LRInfo *OutValue);






%typemap(python,arginit) struct XX * {
	$1 = NULL;
}

%typemap(python,freearg) struct XX * {
	int i;
	if($1->data != NULL) {
	  for (i = 0; i < $1->nn; ++i) {
		free($1->data[i]);
	  }
	  free($1->data);
	}
	free($1);
}

%typemap(python,ignore) struct XX *OutValue {
    $1 = (struct XX *)malloc(sizeof(struct XX));
}

%typemap(python,argout) struct XX *OutValue {
   int i,j;
   PyObject *o, *p, *q;

   o = PyList_New($1->nn);

   for(i = 0; i < $1->nn; ++i) {
	 p = PyList_New($1->k);
	 for(j = 0; j < $1->k; ++j) {
		q = PyFloat_FromDouble($1->data[i][j]);
		PyList_SetItem(p, j, q);
	 }
	 PyList_SetItem(o, i, p);
   }

   if ((!$result) || ($result == Py_None)) {
	   $result = o;
   } else {
	   if (!PyList_Check($result)) {
		   PyObject *o2 = $result;
		   $result = PyList_New(0);
		   PyList_Append($result,o2);
		   Py_XDECREF(o2);
	   }
	   PyList_Append($result,o);
	   Py_XDECREF(o);
   }
}

%typemap(python,in) struct XX * {
 /* Check if is a list */
 if (PyList_Check($input)) {
   int i, j, size, msize;

   size = PyList_Size($input);
   msize = 0;

   $1 = (struct XX*) malloc(sizeof(struct XX));
   $1->data = NULL;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem($input,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o);
		 if ($1->data == NULL) {
			 msize = zsize;
			 $1->nn = size;
			 $1->k = msize;
			 $1->data = (double **)malloc(sizeof(double *)*size);
			 for(j = 0; j < size; ++j) {
				$1->data[j] = (double *)malloc(sizeof(double)*msize);
		  	 }
			 
		 }
		 if (zsize == msize) {
			 for (j = 0; j < zsize; ++j) {
			     PyObject *p = PyList_GetItem(o,j);
				 if (PyFloat_Check(p)) {
					$1->data[i][j] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[i][j] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"sublists must contain doubles");
					 if ($1->data != NULL) {
			 			 for(j = 0; j < size; ++j) {
							free($1->data[j]);
		  				}
		  				free($1->data);
		  			 }
					 free($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"sublists must be of equal size");
			 if ($1->data != NULL) {
				 	for(j = 0; j < size; ++j) {
					free($1->data[j]);
			  	}
			  	free($1->data);
		  	 }
			 free($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"list must contain lists");
		if ($1->data != NULL) {
			for(j = 0; j < size; ++j) {
				free($1->data[j]);
			}
			free($1->data);
	   }
       free($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"not a list");
   return NULL;
 }
}


void Computer(struct XX *input, struct XX *OutValue);














struct NBInfo *NBprepare(struct NBInput *input);
void NBkill(struct NBInfo *p);
void NBcleanup(struct NBInput *p);
void NBquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);
void TANquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);
void NBdivergence(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);
void NBqualityW(struct NBInfo *in, double *w, struct NBModel *m, struct NBResult *OutValue);
void NBsaveScores(struct NBInfo *in);
void NBrememberScores(struct NBInfo *in);
void NBcompareScores(struct NBInfo *in, struct NBResult *OutValue);
void NBexportScores(struct NBInfo *in, int mode, struct NBList *OutValue);
void NBexportProbabilities(struct NBInfo *in, int mode, struct NBList *OutValue);
double NBcompareLists(int n, double *a, double *b);
void NBstoreModel(struct NBInfo *in, double *w, struct NBModel *m);
void NBclassify(struct NBInfo *in, int *ex, struct NBList *OutValue);
void NBclassifyW(struct NBInfo *in, int *ex, struct NBList *OutValue);
void NBupdate(struct NBInfo *in, int attribute, int card, int *values);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

%typemap(python,in) struct KInput * {
 int i, j, size, msize;
 PyObject *arr,*cards;
 /* Check if is a list */

 arr = PyTuple_GetItem($input,0);
 if (PyList_Check(arr)) {

   size = PyList_Size(arr);
   msize = 0;
	//printf("lr:-1\n");
   $1 = (struct KInput *) malloc(sizeof(struct KInput));
   $1->data = NULL;
   $1->card = NULL;
   for (i = 0; i < size; i++) {
     PyObject *o = PyList_GetItem(arr,i);
     if (PyList_Check(o)) {
		 int zsize = PyList_Size(o);
		 if ($1->data == NULL) {
			 msize = zsize;
			 $1->nn = size;
			 $1->k = msize;
			 $1->na = $1->k-1; // number of attributes
			 $1->data = (int **)malloc(sizeof(int *)*(size));
			 for (j = 0; j < size; ++j) {
				$1->data[j] = (int *)malloc(sizeof(int)*($1->k));
			 }
			 $1->card = (int *)malloc(sizeof(int)*($1->k));
		 }
		 if (zsize == msize) {		 				
			 // fetch the attribute values
			 for (j = 0; j < $1->k; ++j) {
			     PyObject *p = PyList_GetItem(o,j);
				 if (PyFloat_Check(p)) {
					// correct +1 because that array starts at 1
					$1->data[i][j] = PyFloat_AsDouble(p);
				 } else if (PyInt_Check(p)) {
					$1->data[i][j] = (double)PyInt_AsLong(p);
				 } else {
					 PyErr_SetString(PyExc_TypeError,"examples must contain doubles or ints");
					 Kcleanup($1);
					 return NULL;
				 }
			 }
		 } else {
			 PyErr_SetString(PyExc_TypeError,"examples must be of equal size");
			 Kcleanup($1);
			 return NULL;
		 }
     } else {
       PyErr_SetString(PyExc_TypeError,"example table must contain examples as lists");
       Kcleanup($1);
       return NULL;
     }
   }
 } else {
   PyErr_SetString(PyExc_TypeError,"tuple element #0 not a list");
   return NULL;
 }


 // fetch the cardinalities
 cards = PyTuple_GetItem($input,1);
 if (!PyList_Check(cards) && PyList_Size(cards) != $1->k) {
   PyErr_SetString(PyExc_TypeError,"tuple element #1 not an attribute cardinality list");
   return NULL;
 }

 for (j = 0; j < $1->k; ++j) {
	 PyObject *p = PyList_GetItem(cards,j);
	 if (PyInt_Check(p)) {
		$1->card[j] = PyInt_AsLong(p);
	 } else {
		 PyErr_SetString(PyExc_TypeError,"card list must contain ints");
		 Kcleanup($1);
		 return NULL;
	 }
 }
}

%typemap(python,ignore) struct KList *OutValue {
    $1 = (struct KList *)malloc(sizeof(struct KList));
}

%typemap(python,argout) struct KList *OutValue {
   int i;
   PyObject *o, *q;

   o = PyList_New($1->n);

   for(i = 0; i < $1->n; ++i) {
	 q = PyFloat_FromDouble($1->l[i]);
	 PyList_SetItem(o, i, q);
   }
   free($1->l);
   free($1);

   if ((!$result) || ($result == Py_None)) {
	   $result = o;
   } else {
	   if (!PyList_Check($result)) {
		   PyObject *o2 = $result;
		   $result = PyList_New(0);
		   PyList_Append($result,o2);
		   Py_XDECREF(o2);
	   }
	   PyList_Append($result,o);
	   Py_XDECREF(o);
   }
}

%typemap(python,ignore) struct KMatrix *OutValue {
    $1 = (struct KMatrix *)malloc(sizeof(struct KMatrix));
}

%typemap(python,argout) struct KMatrix *OutValue {
   int i,j;
   PyObject *o, *p, *q;

   o = PyList_New($1->rows);

   for(i = 0; i < $1->rows; ++i) {
	 p = PyList_New($1->columns);
	 for(j = 0; j < $1->columns; ++j) {
		q = PyFloat_FromDouble($1->l[i*$1->columns + j]);
		PyList_SetItem(p, j, q);
	 }
	 PyList_SetItem(o, i, p);
   }
   free($1->l);
   free($1);

   if ((!$result) || ($result == Py_None)) {
	   $result = o;
   } else {
	   if (!PyList_Check($result)) {
		   PyObject *o2 = $result;
		   $result = PyList_New(0);
		   PyList_Append($result,o2);
		   Py_XDECREF(o2);
	   }
	   PyList_Append($result,o);
	   Py_XDECREF(o);
   }
}


%typemap(python,in) struct KModel * {
	PyObject *p;

	$1 = (struct KModel *) malloc(sizeof(struct KModel));
	$1->nogroups = 0;
	$1->groups = NULL;

	if (PyList_Check($input)) {
		int i, j, msize;

		$1->nogroups = PyList_Size($input);
		if($1->nogroups > 0) {
			$1->groups = (struct KGroup *)malloc(sizeof(struct KGroup)*$1->nogroups);
			for (i = 0; i < $1->nogroups; i++) {
				PyObject *o = PyList_GetItem($input,i);
				PyObject *list = PyTuple_GetItem(o,1);
				msize = PyList_Size(list);
				$1->groups[i].n = msize;
				$1->groups[i].l = (int *)malloc(sizeof(int)*msize);
				for (j = 0; j < msize; ++j)
					$1->groups[i].l[j] = PyInt_AsLong(PyList_GetItem(list,j));
				$1->groups[i].times = PyInt_AsLong(PyTuple_GetItem(o,0));
			}
		}
	} else {
		PyErr_SetString(PyExc_TypeError,"tuple element #1 not a list of groups");
		return NULL;
	}
}

%typemap(python,arginit) struct KModel * {
	$1 = NULL;
}

%typemap(python,freearg) struct KModel * {
	int i;
	if($1->groups != NULL) {
	  for (i = 0; i < $1->nogroups; ++i) {
		free($1->groups[i].l);
	  }
	  free($1->groups);
	}
	free($1);
}

/////

%typemap(python,in) struct KModels * {
	PyObject *p;

	$1 = (struct KModels *) malloc(sizeof(struct KModels));
	$1->nomodels = 0;
	$1->models = NULL;

	if (PyList_Check($input)) {
		int i, j, k,msize;
		KModel *mod;

		$1->nomodels = PyList_Size($input);
		if($1->nomodels > 0) {
			$1->models = (struct KModel *)malloc(sizeof(struct KModel)*$1->nomodels);
			for (k = 0; k < $1->nomodels; k++) {
				p = PyList_GetItem($input,k);
				mod = $1->models + k;
				mod->groups = NULL;
				mod->nogroups = PyList_Size(p);
				if(mod->nogroups > 0) {
					mod->groups = (struct KGroup *)malloc(sizeof(struct KGroup)*mod->nogroups);
					for (i = 0; i < mod->nogroups; i++) {
						PyObject *o = PyList_GetItem(p,i);
						PyObject *list = PyTuple_GetItem(o,1);
						msize = PyList_Size(list);
						mod->groups[i].n = msize;
						mod->groups[i].l = (int *)malloc(sizeof(int)*msize);
						for (j = 0; j < msize; ++j)
							mod->groups[i].l[j] = PyInt_AsLong(PyList_GetItem(list,j));
						mod->groups[i].times = PyInt_AsLong(PyTuple_GetItem(o,0));
					}
				}
			}
		}
	} else {
		PyErr_SetString(PyExc_TypeError,"tuple element #1 not a list of groups");
		return NULL;
	}
}

%typemap(python,arginit) struct KModels * {
	$1 = NULL;
}

%typemap(python,freearg) struct KModels * {
	int i,j,k;
	KModel *mod;
	if($1->models != NULL) {
		for (k = 0; k < $1->nomodels; ++k) {
			mod = $1->models+k;
			if(mod->groups != NULL) {
				for (i = 0; i < mod->nogroups; ++i) {
					free(mod->groups[i].l);
				}
				free(mod->groups);
			}
		}
		free($1->models);
	}
	free($1);
}


%typemap(python,in) struct KArray * {
	PyObject *p;

	$1 = (struct KArray *) malloc(sizeof(struct KArray));
	$1->n = 0;
	$1->l = NULL;

	if (PyList_Check($input)) {
		int i, j, msize;

		$1->n = PyList_Size($input);
		if($1->n > 0) {
			$1->l = (double *)malloc(sizeof(double)*$1->n);
			for (i = 0; i < $1->n; i++) {
				$1->l[i] = PyFloat_AsDouble(PyList_GetItem($input,i));
			}
		}
	} else {
		PyErr_SetString(PyExc_TypeError,"this is not a list of doubles");
		return NULL;
	}
}

%typemap(python,arginit) struct KArray * {
	$1 = NULL;
}

%typemap(python,freearg) struct KArray * {
	int i;
	if($1->l != NULL) {
	  free($1->l);
	}
	free($1);
}

void Ksetmodel(struct KInfo *in, struct KModel *m);
void Kaddmodel(struct KInfo *in, struct KModel *m);
void Ktestaddition(struct KInfo *in, struct KModel *m, struct KList *OutValue);

void Kdie(struct KInfo *p);
void Kuse(struct KInfo *in, int *ex, struct KList *OutValue);
double Kvalidate(struct KInfo *in, int *ex);
struct KInfo *Kremember(struct KInput *input, double prior);
void Klearn(struct KInfo *in, int samples, int depth, struct KMatrix *OutValue);
void KgetDOF(struct KInfo *in, struct KList *OutValue);
void Kcheckreversal(struct KInfo *in, struct KList *OutValue);
void Ksetensemble(struct KInfo *in, struct KModels *ms, struct KArray *weights);
void Kuseensemble(struct KInfo *in, int *ex, struct KList *OutValue);
void Ktestmodels(struct KInfo *in, struct KModels *ms, int samples, struct KMatrix *OutValue);
