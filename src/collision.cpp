#include <wmrde/collision.h>
#include <wmrde/algebra/rotation.h>
#include <wmrde/util/rootfinding.h>
#include <wmrde/util/rosout.h>

namespace wmrde
{

//return transform from contact to world coords
HTransform calcHTContactToWorld(
    const Vec3& point, //contact point (in world coords)
    const Vec3& normal, //ground surface normal (in world coords)
    const HTransform& HT_wheel_to_world) //transform from wheel to world coords
{
  //cross the wheel frame y-axis with the normal to get
  //contact frame x-axis that is orthogonal to both
  Vec3 x_axis = (HT_wheel_to_world.R.col(1).cross(normal)).normalized();

  HTransform out;
  out.R.col(0) = x_axis;
  out.R.col(1) = normal.cross(x_axis);
  out.R.col(2) = normal;
  out.t = point;

  return out;
}

void calcContactFrameDiscretize(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactFrame& contact)
{

  int tangent_idx = -1;
  Real dz = RealInf();
  Vec3 point;
  Vec3 normal;

  for (size_t i = 0; i < wheel_geom.tangents.size(); i++)
  {
    Vec3 this_point = HT_wheel_to_world.applyTo(wheel_geom.tangents[i].point);
    Real this_dz;
    Vec3 this_normal;
    surfs.getDistance(this_point, this_dz, this_normal);
    if (this_dz < dz)
    {
      tangent_idx = i;
      dz = this_dz;
      point = this_point;
      normal = this_normal;
    }
  }

  if (tangent_idx >= 0)
  {
    contact.dz = dz;
    contact.wheel_tangent_dir = wheel_geom.tangents[tangent_idx].direction;
    contact.HT_world = calcHTContactToWorld(point, normal, HT_wheel_to_world);
    contact.HT_wheel = HT_wheel_to_world.composeInvWith(contact.HT_world);
  }
}

void calcContactFrameRootfind(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactFrame& contact)
{
  //find contact angle where dot product of tangent to wheel and surface normal == 0

  //get range of contact angles
  Real angle_range[2];
  {
    Real xproj = HT_wheel_to_world.R(2,0); //dot(world z axis, wheel x axis)
    Real zproj = HT_wheel_to_world.R(2,2); //dot(world z axis, wheel z axis)
    Real mid = atan2(xproj,zproj); // middle of range

    Real half_range = M_PI/2.0;

    angle_range[0] = mid - half_range;
    angle_range[1] = mid + half_range;
  }
  ROS_DEBUG("angle range [%f, %f]", angle_range[0], angle_range[1]); //DEBUGGING

  Real tolx = degToRad(.5);
  Real tolfx = 1e-2;

  //set in lambda
  Real dz;
  Vec3 point; //contact point in world coords
  Vec3 normal; //surface normal in world coords
  Tangent tangent;

  //lambda closure, requires C++11
  auto dotProductTangentNormal = [&] (const Real angle) -> Real
  {
    //convert contact angle to tangent
    tangent = contactAngleToTangent(wheel_geom.radius, angle); //in wheel coords

    //convert tangent to world coords
    point = HT_wheel_to_world.applyTo(tangent.point);
    Vec3 tangent_dir = HT_wheel_to_world.R*tangent.direction; //rotate into world coords

    //get surface normal
    surfs.getDistance(point, dz, normal);

    Real dp = tangent_dir.dot(normal);
    ROS_DEBUG("  angle = %f, tangent.dot(normal) = %f", angle, dp);

    return dp;
  };

  //uncomment *one* of the following:
  Real angle;
//  angle = findRootBisection(angle_range[0],angle_range[1],dotProductTangentNormal,tolx,tolfx);
  angle = findRootBrents(angle_range[0],angle_range[1],dotProductTangentNormal,tolx,tolfx); //faster

  contact.dz = dz;
  contact.wheel_tangent_dir = tangent.direction;
  contact.HT_world = calcHTContactToWorld(point, normal, HT_wheel_to_world);
  contact.HT_wheel = HT_wheel_to_world.composeInvWith(contact.HT_world);
}

} //namespace
