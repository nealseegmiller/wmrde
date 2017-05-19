#ifndef _WMRDE_COLLISION_H_
#define _WMRDE_COLLISION_H_

#include <wmrde/surface/surface.h>
#include <wmrde/contactgeom.h>

namespace wmrde
{

/*!
 * Convert contact angle to point. Convention: angle 0 corresponds
 * to point [0,0,-radius] and positive angle is ccw about y axis.
 * \param radius The wheel radius
 * \param angle The contact angle
 * \return The contact point in wheel coords
 */
inline Vec3 contactAngleToPoint(const Real radius, const Real angle)
{
  Vec3 pt;
  pt << -sin(angle)*radius, 0.0, -cos(angle)*radius;
  return pt;
}

void setWheelGeomPoints(
    const double min_angle,
    const double max_angle,
    const int num_points,
    WheelGeom& wheel_geom) //angle must be set
{
  Real angle_step = (max_angle - min_angle)/(num_points-1);
  wheel_geom.points.resize(num_points);
  for (int i = 0; i < num_points; i++)
  {
    Real angle = min_angle + i*angle_step;
    wheel_geom.points[i] = contactAngleToPoint(wheel_geom.radius, angle);
  }
}

/*!
 * Determine the contact point between wheel and surface geometries
 * using discretization.
 * \param surfs The surfaces
 * \param wheel_goem The wheel geometry, requires that points be set
 * \param HT_wheel_to_world The transform from wheel to world coords
 * \param contact The contact geometry output.
 * \return index of wheel_geom point that was selected as contact point
*/
int calcContactGeomDiscretization(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactGeom& contact);

/*!
 * Determine the contact point between wheel and surface geometries
 * using rootfinding.
 * \param surfs The surfaces
 * \param wheel_goem The wheel geometry, requires radius be set
 * \param HT_wheel_to_world The transform from wheel to world coords
 * \param contact The contact geometry output.
*/
void calcContactGeomRootfinding(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactGeom& contact);

} //namespace

#endif //_WMRDE_COLLISION_H_
