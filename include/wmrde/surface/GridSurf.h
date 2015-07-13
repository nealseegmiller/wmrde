//GridSurf.h
//The GridSurf class is derived from the Surface base class.
//This class specifies a uniformly-spaced grid of height data. Interpolation is required to compute heights and surface normals.

#ifndef _WMRSIM_GRIDSURF_H_
#define _WMRSIM_GRIDSURF_H_

//to read in data from csv file
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

#include <external/xs_Float.h> //cast float to int

#include <surface/Surface.h>
#include <algebra/transform.h>
#include <common/interp.h>

//need public so conversion from GridSurf* to Surface* is possible?
class GridSurf : public Surface {
private:
	static const int buffer = 201*201; //size of Zdata

	//x,y
	Real lowerlim[2]; //lower limit
	Real upperlim[2]; //upper limit

	int n[2]; //number

	Real Z[buffer]; //height data

	//DEPENDENT
	Real d[2]; //spacing
	//precompute gradients to speed up surface normal calculation
	Real DZDX[buffer];
	Real DZDY[buffer];

	void calcGradients();
	int get_loc(const Vec3 pt);

public:

	const Real get_lowerlimx() const { return lowerlim[0]; }
	const Real get_lowerlimy() const { return lowerlim[1]; }
	const int get_nx() const { return n[0]; }
	const int get_ny() const { return n[1]; }
	const Real get_dx() const { return d[0]; }
	const Real get_dy() const { return d[1]; }
	Real* get_Z() { return Z; }

	GridSurf() {} //constructor
	int readfile(const std::string filename);

	//need to re-declare virtual functions?
	//compiler doesn't check if signature is the same?
	int surfaceHeight(const Vec3 pt, int loc, Real& height);
	int surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal);
	int surfaceNormal(const Vec3 pt, int loc, Vec3 normal);
	
};


#endif //_WMRSIM_GRIDSURF_H_