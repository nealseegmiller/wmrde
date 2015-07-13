//interp.h
//functions for linear, bilinear, and trilinear interpolation

#ifndef _WMRDE_INTERP_H_
#define _WMRDE_INTERP_H_

#include <wmrde/util/common_util.h>
#include <wmrde/util/xs_Float.h> //cast float to int

#define S2I2(i,j,n1) (i + (j)*(n1))
#define S2I3(i,j,k,n1,n2) (i + (j)*(n1) + (k)*(n1)*(n2))

inline void Interp_xn ( const int dim, const Real lim[], const Real d[], const Real x[], int I[], Real xn[] ) {
  //dim:    dimension 1,2,3
  //lim:    lower limit for x
  //d:    spacing
  //x:
  //I:    index
  //xn:   x normalized between 0 and 1

  for (int i=0; i<dim; i++) {
    //I[i] = (int) floor((x[i]-lim[i])/d[i]);
    I[i] = xs_FloorToInt((x[i]-lim[i])/d[i]); //faster?
    xn[i] = (x[i] - (lim[i]+I[i]*d[i])) / d[i]; //rescale between 0 and 1
  }
}

inline void linearInterp_xn( const Real lim[], const Real d[], const Real x[], int I[], Real xn[] ) {
  Interp_xn( 1, lim, d, x, I, xn );
}
inline void bilinearInterp_xn( const Real lim[], const Real d[], const Real x[], int I[], Real xn[] ) {
  Interp_xn( 2, lim, d, x, I, xn );
}
inline void trilinearInterp_xn( const Real lim[], const Real d[], const Real x[], int I[], Real xn[] ) {
  Interp_xn( 3, lim, d, x, I, xn );
}

inline void linearInterp_f(const Real* F, const int I[1], Real f[2]) {
  f[0] = F[I[0]];
  f[1] = F[I[0]+1];
}
inline void bilinearInterp_f(const Real* F, const int I[2], const int n[1], Real f[4]) {
  f[0] = F[S2I2(I[0],   I[1], n[0])];
  f[1] = F[S2I2(I[0]+1, I[1], n[0])];
  f[2] = F[S2I2(I[0],   I[1]+1, n[0])];
  f[3] = F[S2I2(I[0]+1, I[1]+1, n[0])];
}
inline void trilinearInterp_f(const Real* F, const int I[3], const int n[2], Real f[8]) {
  f[0] = F[S2I3(I[0],   I[1], I[2], n[0], n[1])];
  f[1] = F[S2I3(I[0]+1, I[1], I[2], n[0], n[1])];
  f[2] = F[S2I3(I[0],   I[1]+1, I[2], n[0], n[1])];
  f[3] = F[S2I3(I[0]+1, I[1]+1, I[2], n[0], n[1])];
  f[4] = F[S2I3(I[0],   I[1], I[2]+1, n[0], n[1])];
  f[5] = F[S2I3(I[0]+1, I[1], I[2]+1, n[0], n[1])];
  f[6] = F[S2I3(I[0],   I[1]+1, I[2]+1, n[0], n[1])];
  f[7] = F[S2I3(I[0]+1, I[1]+1, I[2]+1, n[0], n[1])];
}


inline Real linearInterp( const Real f[2], const Real xn[1]) {
  return f[0] + xn[0]*(f[1]-f[0]);
}

inline Real bilinearInterp( const Real f[4], const Real xn[2]) {
  Real f0 = linearInterp( f+0, xn+0 );
  Real f1 = linearInterp( f+2, xn+0 );
  return f0 + xn[1]*(f1-f0);
}

inline Real trilinearInterp( const Real f[8], const Real xn[3]) {
  Real f0 = bilinearInterp ( f+0, xn+0 );
  Real f1 = bilinearInterp ( f+4, xn+0 );
  return f0 + xn[2]*(f1-f0);
}

#endif //_WMRDE_INTERP_H_
