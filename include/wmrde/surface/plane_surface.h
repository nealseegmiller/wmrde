#ifndef _WMRDE_PLANE_SURFACE_H_
#define _WMRDE_PLANE_SURFACE_H_

#include <wmrde/surface/surface.h>
#include <wmrde/algebra/transform.h>

namespace wmrde
{

struct BoundingBox
{
  Real x_min;
  Real x_max;
  Real y_min;
  Real y_max;

  BoundingBox() {}
  BoundingBox(
      Real _x_min,
      Real _x_max,
      Real _y_min,
      Real _y_max)
    :
      x_min(_x_min),
      x_max(_x_max),
      y_min(_y_min),
      y_max(_y_max)
  {}

  bool in(Real x, Real y) const
  {
    return x > x_min && x < x_max &&
           y > y_min && y < y_max;
  }
  bool inOrOn(Real x, Real y) const
  {
    return x >= x_min && x <= x_max &&
           y >= y_min && y <= y_max;
  }
};

class PlaneSurface : public Surface
{
 public:
	//overload constructor
  PlaneSurface() : initialized_(false), bounded_(false) {}
  PlaneSurface(const HTransform& T) : bounded_(false) { setTransform(T); }
  PlaneSurface(const HTransform& T, const BoundingBox& bb)
  {
    setTransform(T);
    setBoundingBox(bb);
  }

  bool initialized() const { return initialized_; }
  bool bounded() const { return bounded_; }

  /*!
   * Set the plane transform.
   * \param T The transform from plane to world coordinates: T.R.col(2) is the plane normal vector and T.t is the plane origin
   */
  void setTransform(const HTransform& T)
  {
    T_ = T;
    Tinv_ = T.inverse(); //transform from world to plane coordinates
    plane_eq = { Tinv_.R(2,0), Tinv_.R(2,1), Tinv_.R(2,2), Tinv_.t(2) }; //bottom row of Tinv
    initialized_ = true;
  }
  /*!
   * Set the bounding box. for bounded planes, the getHeight, etc. functions will
   * return false if the query point is outside of the bounding box
   * \param bb The bounding box in plane coordinates
   */
  void setBoundingBox(const BoundingBox& bb)
  {
    bb_ = bb;
    bounded_ = true;
  }

  //override virtual functions in Surface class
	bool getHeight(const Vec3& pt, Real& height) const override;
	bool getNormal(const Vec3& pt, Vec3& normal) const override;
	bool getDistance(const Vec3& pt, Real& distance, Vec3& normal) const override;
	
 private:
	bool initialized_;
	bool bounded_;
	BoundingBox bb_;
	HTransform T_; //transform from plane to world coordinates
	HTransform Tinv_; //transform from world to plane coordinates
	std::array<Real,4> plane_eq; //plane equation coefficients: a*x + b*y + c*z + d = 0

	//return true if point is inOrOn bounding box
	bool inBounds(const Vec3& pt) const ;

	friend class WmrdeRosInterface;
};

} //namespace

#endif //_WMRDE_PLANE_SURFACE_H_
