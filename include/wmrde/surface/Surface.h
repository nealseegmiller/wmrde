//Surface.h
//abstract Surface class
//all pure (=0) virtual functions must be defined by derived classes

#ifndef _WMRDE_SURFACE_H_
#define _WMRDE_SURFACE_H_

namespace wmrde
{

} //namespace

#endif //_WMRDE_SURFACE_H_

/*
#include <vector>
#include <memory>
#include <wmrde/algebra/linalg3.h>

class Surface {
public:
	//pt:			x,y location of points (z not required)
	//loc:			location (triangle index for TriMeshSurf). if < 0 find loc
	//height:		surface height at point, NaN if out of bounds
	//(return)		loc, < 0 if out of bounds
	virtual int surfaceHeight(const Vec3 pt, int loc, Real& height) = 0;

	//pt:			x,y,z location of point
	//dz:			distance of point from the surface (measured along the surface normal). negative if below surface, NaN if out of bounds
	//normal:		surface normal, may be null
	virtual int surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal) = 0;

	//pt:			x,y location of points (z not required)
	//normal:		surface normal (unit vector) at point
	virtual int surfaceNormal(const Vec3 pt, int loc, Vec3 normal) = 0;
};

//a vector of pointers to dynamically allocated Surface objects
//use unique_ptr (C++11) for simple memory managment
typedef std::vector<std::unique_ptr<Surface>> SurfaceVector;

//todo, should also output loc?

//& to pass by reference
//return surface index (in vector)
int surfacesHeight(const SurfaceVector& surfaces, const Vec3 pt, Real& height);
int surfacesDz(const SurfaceVector& surfaces, const Vec3 pt, Real& dz, Vec3 normal);
*/


