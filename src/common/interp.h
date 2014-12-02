//interp.h
//functions for linear, bilinear, and trilinear interpolation

#ifndef _WMRSIM_INTERP_H_
#define _WMRSIM_INTERP_H_

#include <common/common.h>
#include <external/xs_Float.h> //cast float to int

#define S2I2(i,j,n1) (i + (j)*(n1))
#define S2I3(i,j,k,n1,n2) (i + (j)*(n1) + (k)*(n1)*(n2))

inline void Interp_xn ( const int dim, const Real lim[], const Real d[], const Real x[], int I[], Real xn[] );

void linearInterp_xn( const Real lim[], const Real d[], const Real x[], int I[], Real xn[] );
void bilinearInterp_xn( const Real lim[], const Real d[], const Real x[], int I[], Real xn[] );
void trilinearInterp_xn( const Real lim[], const Real d[], const Real x[], int I[], Real xn[] );

void linearInterp_f(const Real* F, const int I[1], Real f[2]);
void bilinearInterp_f(const Real* F, const int I[2], const int n[1], Real f[4]);
void trilinearInterp_f(const Real* F, const int I[3], const int n[2], Real f[8]);

Real linearInterp( const Real f[2], const Real xn[1] );
Real bilinearInterp( const Real f[4], const Real xn[2] );
Real trilinearInterp( const Real f[8], const Real xn[3] );


#endif