%module Skin_Setup
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
#include "Skin_Setup.h"
%}

%include "numpy.i"
%init %{
  import_array()
%}

%include "cpointer.i"
%pointer_functions(double, doubleP)

%include "Skin.h"
%include "Skin_Setup.h"

