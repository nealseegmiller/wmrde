//common.h
//#defines used throughout the project

#ifndef _WMRDE_COMMON_UTIL_H_
#define _WMRDE_COMMON_UTIL_H_

//#include <iostream>
#include <cmath>
#include <limits> //for NaN
#include <algorithm> //for std::max()

#include <wmrde/common.h>
#include <wmrde/options.h> //for typedef Real


//TODO, move this to options.h?
//defines based on options

//now defined in common.h
//#if WMRSIM_DOUBLE_PRECISION
//typedef double Real;
//#else
//typedef float Real;
//#endif

#if WMRSIM_USE_QUATERNION
#define SIZEORIENT 4
#else
#define SIZEORIENT 3
#endif

//for output only in Debug mode
#ifdef _DEBUG
#define DEBUG_CERR(x) std::cerr << x << std::endl;
#else
#define DEBUG_CERR(x) //do {} while (0)
#endif

//signum
#define TYPESIGN(val,Type) ( Type((val > 0) - (val < 0)) )
#define REALSIGN(val) ( Real((val > 0) - (val < 0)) )

//for Real
#define REALRAND ( Real(rand()) / RAND_MAX ) //pseudo-random number between 0 and 1, Real precision
#define REALNAN std::numeric_limits<Real>::quiet_NaN()
#define REALMAX std::numeric_limits<Real>::max()

//for angles
#undef M_PI
#define M_PI Real(3.141592653589793)

//#define M_TAU (2*M_PI)
const Real M_TAU = 2*M_PI; //faster?

//convert angles between degrees and radians
#define DEGTORAD(a) ((a)*M_PI/180)
#define RADTODEG(a) ((a)*180/M_PI)
//wrap angle to range of 2*pi radians, (360 degrees)
//c++ fmod() is not equivalent to Matlab mod(), round towards zero vs. floor
#define WRAPRAD(a,max) ( a+max - floor((a+max)/M_TAU)*M_TAU + (max-M_TAU) )
#define WRAPDEG(a,max) ( a+max - floor((a+max)/360)*360 + (max-360) )
//difference two angles, a-b, range of -pi to pi radians (-180 to 180 degrees)
#define DIFFRAD(a,b) ( WRAPRAD((a-b),M_PI) )
#define DIFFDEG(a,b) ( WRAPDEG((a-b),180) )


#endif  //_WMRDE_COMMON_H_
