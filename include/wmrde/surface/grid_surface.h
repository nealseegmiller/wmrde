//GridSurf.h
//The GridSurf class is derived from the Surface base class.
//This class specifies a uniformly-spaced grid of height data. Interpolation is required to compute heights and surface normals.

#ifndef _WMRDE_GRID_SURFACE_H_
#define _WMRDE_GRID_SURFACE_H_

#include <string>
#include <assert.h>

#include <wmrde/surface/surface.h>
#include <wmrde/util/gridded_interpolant.h>

namespace wmrde
{

//need public so conversion from GridSurf* to Surface* is possible?
class GridSurface : public Surface {
 private:

//	GriddedInterpolant<Real,2> Gh_; //TODO, for faster height only calculation, but uses more memory
	GriddedInterpolant<Vec3,2> G_; //for height and surface normal

 public:

	GridSurface() {}
	GridSurface(
	    const GridVector& X,
	    const GridVector& Y,
	    const std::vector<Real>& Z); //elevations

	GridSurface(const std::string& filename); //TODO

  //override virtual functions in Surface class
  bool getHeight(const Vec3& pt, Real& height) const override;
  bool getNormal(const Vec3& pt, Vec3& normal) const override;
  bool getDistance(const Vec3& pt, Real& distance, Vec3& normal) const override;
	
};

} //namespace

#endif
