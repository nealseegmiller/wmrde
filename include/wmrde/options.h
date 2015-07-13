//options.h
//compile-time options for WMR simulations

#ifndef _WMRDE_OPTIONS_H_
#define _WMRDE_OPTIONS_H_

#include <string>

//OPTIONS, set these flags to 1 or 0
#define WMRSIM_USE_QUATERNION 0 //else use Euler angles
#define WMRSIM_DOUBLE_PRECISION 1 //else single (float) precision
//single precision is slower?! lots of "trunction, possible loss of data" warnings
#define WMRSIM_ENABLE_ANIMATION 1 //include WmrAnimation, OGRE dependencies

inline std::string ResourceDir() {
	//return std::string("C:/Users/nseegmil.NREC-014635/Documents/Thesis/src/resource/");
  return std::string("/home/neal/Projects/wmrde/resource/");
}

inline std::string CADdir() {
	//return std::string("C:/Users/nseegmil.NREC-014635/Dropbox/CAD/");
  return std::string("/home/neal/Projects/wmrde/CAD/");
}


//don't modify below here

/*
//defines based on options
#if WMRSIM_DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

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
*/

#endif  //_WMRDE_OPTIONS_H_
