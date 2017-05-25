#ifndef _WMRDE_COLLISION_H_
#define _WMRDE_COLLISION_H_

#include <wmrde/surface/surface.h>
#include <wmrde/contactframe.h>
#include <wmrde/wheelgeom.h>

namespace wmrde
{

/*!
 * Convert contact angle to point. Convention: angle 0 corresponds
 * to point [0,0,-radius] and positive angle is ccw about y axis.
 * \param radius The wheel radius
 * \param angle The contact angle
 * \return The contact point in wheel coords
 */
inline Tangent contactAngleToTangent(const Real radius, const Real angle)
{
  Tangent out;
  out.point << -sin(angle)*radius, 0.0, -cos(angle)*radius;
  out.direction << -cos(angle), 0.0, sin(angle);
  return out;
}

inline void discretizeWheelGeom(
    const double min_angle,
    const double max_angle,
    const int num_points,
    WheelGeom& wheel_geom) //angle must be set
{
  Real angle_step = (max_angle - min_angle)/(num_points-1);
  wheel_geom.tangents.resize(num_points);
  for (int i = 0; i < num_points; i++)
  {
    Real angle = min_angle + i*angle_step;
    wheel_geom.tangents[i] = contactAngleToTangent(wheel_geom.radius, angle);
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
void calcContactFrameDiscretize(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactFrame& contact);

/*!
 * Determine the contact point between wheel and surface geometries
 * using rootfinding.
 * \param surfs The surfaces
 * \param wheel_goem The wheel geometry, requires radius be set
 * \param HT_wheel_to_world The transform from wheel to world coords
 * \param contact The contact geometry output.
*/
void calcContactFrameRootfind(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactFrame& contact);

} //namespace

#endif //_WMRDE_COLLISION_H_
