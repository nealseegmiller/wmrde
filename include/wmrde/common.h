//common.h

#ifndef _WMRDE_COMMON_H_
#define _WMRDE_COMMON_H_

#include <limits>

//choose double or single precision
typedef double Real;
//typedef float Real;

inline Real RealNan() { return std::numeric_limits<Real>::quiet_NaN(); }
inline Real RealInf() { return std::numeric_limits<Real>::infinity(); }

#endif
