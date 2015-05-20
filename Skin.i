%module Skin
%{
#define SWIG_FILE_WITH_INIT
#include <cvode/cvode.h>
#include <cvode/cvode_band.h>        /* prototype for CVBand */
#include <cvode/cvode_dense.h>   
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_band.h>  /* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS and EXP */
#include "Skin.h"
%}

%include "numpy.i"
%init %{
  import_array()
%}

%apply (double* IN_ARRAY1, int DIM1) {(double *ret, int dim_ret)}
%apply (double* IN_ARRAY1, int DIM1) {(double *coord_x, int dim)};
%apply (double* IN_ARRAY1, int DIM1) {(double *coord_y, int dim)};

%include "Skin.h"
