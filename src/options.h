//options.h
//compile-time options for WMR simulations

#ifndef _WMRSIM_OPTIONS_H_
#define _WMRSIM_OPTIONS_H_

#include <string>

//OPTIONS, set these flags to 1 or 0
#define WMRSIM_USE_QUATERNION 0 //else use Euler angles
#define WMRSIM_DOUBLE_PRECISION 1 //else single (float) precision
//single precision is slower?! lots of "trunction, possible loss of data" warnings
#define WMRSIM_ENABLE_ANIMATION 1 //include WmrAnimation, OGRE dependencies

inline std::string ResourceDir() {
	//return std::string("C:/Users/nseegmil.NREC-014635/Documents/Thesis/src/resource/");
  return std::string("/home/neal/Projects/wmrde/src/resource/");
}

inline std::string CADdir() {
	//return std::string("C:/Users/nseegmil.NREC-014635/Dropbox/CAD/");
  return std::string("/home/neal/Projects/wmrde/CAD/");
}


#endif  //_WMRSIM_OPTIONS_H_
