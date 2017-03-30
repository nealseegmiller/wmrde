#include <wmrde/surface/plane_surface.h>

namespace wmrde
{

bool PlaneSurface::inBounds(const Vec3& pt) const
{
  if (bounded_)
  {
    Vec3 pt_plane = Tinv_.applyTo(pt); //transform point into plane coordinates
    if (!bb_.inOrOn(pt_plane[0], pt_plane[1]))
    {
      return false;
    }
  }
  return true;
}

bool PlaneSurface::getHeight(const Vec3& pt, Real& height) const
{
  if (!initialized_) { return false; }
  height = -(plane_eq[0]*pt[0] + plane_eq[1]*pt[1] + plane_eq[3])/plane_eq[2]; // z = -(a*x + b*y + d)/c
  return inBounds(pt);
}

bool PlaneSurface::getNormal(const Vec3& pt, Vec3& normal) const
{
  if (!initialized_) { return false; }
  for (size_t i = 0; i < 2; i++) { normal[i] = T_.R(i,2); }
  return inBounds(pt);
}

bool PlaneSurface::getDistance(const Vec3& pt, Real& distance, Vec3& normal) const
{
  if (!initialized_) { return false; }

  distance = plane_eq[0]*pt[0] + plane_eq[1]*pt[1] + plane_eq[2]*pt[2] + plane_eq[3];
  for (size_t i = 0; i < 2; i++) { normal[i] = T_.R(i,2); }
	return inBounds(pt);
}


} //namespace

