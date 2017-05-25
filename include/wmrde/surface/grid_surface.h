//GridSurf.h
//The GridSurf class is derived from the Surface base class.
//This class specifies a uniformly-spaced grid of height data. Interpolation is required to compute heights and surface normals.

#ifndef _WMRDE_GRID_SURFACE_H_
#define _WMRDE_GRID_SURFACE_H_

#include <wmrde/surface/surface.h>
#include <wmrde/util/gridded_interpolant.h>

namespace wmrde
{

//need public so conversion from GridSurf* to Surface* is possible?
class GridSurface : public Surface {
 public:
	GridSurface();
	GridSurface(
	    const GridVector& X,
	    const GridVector& Y,
	    const std::vector<Real>& Z); //elevations

	GridSurface(const std::string& filename); //TODO

	/*!
	 * Get elevation grid point by index. Does not interpolate.
	 * \param i the index in X
	 * \param j the index in Y
	 * \return the x,y,z point at the index
	 */
	Vec3 getPoint(const int i, const int j) const ;
	bool initialized() const { return initialized_; }

  //override virtual functions in Surface class
  bool getHeight(const Vec3& pt, Real& height) const override;
  bool getNormal(const Vec3& pt, Vec3& normal) const override;
  bool getDistance(const Vec3& pt, Real& distance, Vec3& normal) const override;
	
 private:
//  GriddedInterpolant<Real,2> Gh_; //TODO, for faster height only calculation, but uses more memory
  GriddedInterpolant<Vec3,2> G_; //for height and surface normal
  bool initialized_;

  friend class WmrdeRosInterface;
};

} //namespace

#endif
