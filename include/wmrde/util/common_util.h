//common.h
//typedef of Real and common mathematical functions
//e.g. for converting angles

#ifndef _WMRDE_COMMON_H_
#define _WMRDE_COMMON_H_

//#include <iostream>
#include <cmath>
#include <limits> //for NaN
#include <algorithm> //for std::max()

#include <wmrde/options.h> //for typedef Real

//
#if WMRDE_USE_DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

//for Real
#define REALRAND ( static_cast<Real>(rand()) / RAND_MAX ) //pseudo-random number between 0 and 1, Real precision
#define REALNAN std::numeric_limits<Real>::quiet_NaN()
#define REALMAX std::numeric_limits<Real>::max()

//signum
inline Real signReal(const Real val) { return (val > 0) - (val < 0); }

//#define M_TAU (2.0*M_PI)
const Real M_TAU = 2.0*M_PI; //faster?

//convert angles between degrees and radians
inline Real degToRad(const double val) { return val*M_PI/180.0; }
inline Real radToDeg(const double val) { return val*180.0/M_PI; }

//wrap angle to range of 2*pi radians,
//e.g. if max is pi, returns value in range -pi to pi
inline Real mod2Pi(const double angle, const double max)
{
  return angle+max - floor((angle+max)/M_TAU)*M_TAU + (max-M_TAU);
}

#endif  //_WMRDE_COMMON_H_
