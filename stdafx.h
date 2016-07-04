// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_band.h>        /* prototype for CVBand */
#include <cvode/cvode_dense.h>   
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_band.h>  /* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS and EXP */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// TODO: reference additional headers your program requires here
#include "except.h"
#define NTHREADS 1 // number of threads for parallel computing
//#define _DEBUG_3_
