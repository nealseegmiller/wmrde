//PlaneSurf.h
//The PlaneSurf class is derived from the surface class
//This class specifies an infinite planar surface.

#ifndef _WMRSIM_PLANESURF_H_
#define _WMRSIM_PLANESURF_H_

#include <assert.h>

#include <surface/Surface.h>
#include <algebra/transform.h>

//need public so conversion from PlaneSurf* to Surface* is possible?
class PlaneSurf : public Surface {
public:
	static const int MAXNV = 10; //max number of vertices
private:
	Real pec[4]; //plane equation coefficients [a b c d]*[x y z 1]^T = 0
	
public:
	//overload constructor
	PlaneSurf(const Real PlaneEqnCoeff[4]);
	PlaneSurf(const Vec3 normal, const Vec3 point);

	//need to re-declare virtual functions?
	//compiler doesn't check if signature is the same?
	int surfaceHeight(const Vec3 pt, int loc, Real& height);
	int surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal);
	int surfaceNormal(const Vec3 pt, int loc, Vec3 normal);
	
};

#endif //_WMRSIM_PLANESURF_H_