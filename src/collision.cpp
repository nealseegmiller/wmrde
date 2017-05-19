#include <wmrde/collision.h>
#include <wmrde/algebra/rotation.h>
#include <wmrde/util/rootfinding.h>
#include <wmrde/rosout.h>

namespace wmrde
{

//calculate transform from contact to world coords.
void calcHTContactToWorld(
    const Vec3& point, //contact point (all inputs in world coords)
    const Vec3& normal, //surface normal
    const Vec3& wheel_yaxis, //wheel frame y-axis
    HTransform& HT_contact_to_world) //output, transform from contact to world coords
{
  Vec3 xaxis = (wheel_yaxis.cross(normal)).normalized();
  Vec3 yaxis = normal.cross(xaxis);

  //copy to HT
  HT_contact_to_world.R.col(0) = xaxis;
  HT_contact_to_world.R.col(1) = yaxis;
  HT_contact_to_world.R.col(2) = normal;
  HT_contact_to_world.t = point;
}

int calcContactGeomDiscretization(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactGeom& contact)
{
  Real min_dz = RealInf();
  int contact_point_idx = -1;
  for (size_t i = 0; i < wheel_geom.points.size(); i++)
  {
    Vec3 point = HT_wheel_to_world.applyTo(wheel_geom.points[i]);
    Real dz;
    Vec3 normal;
    surfs.getDistance(point, dz, normal);
    if (dz < min_dz)
    {
      calcHTContactToWorld(point, normal, HT_wheel_to_world.R.col(1), contact.HT_world);
      min_dz = dz;
      contact_point_idx = i;
    }
  }
  return contact_point_idx;
}

void calcContactGeomRootfinding(
    const Surfaces& surfs,
    const WheelGeom& wheel_geom,
    const HTransform& HT_wheel_to_world,
    ContactGeom& contact)
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
  ROS_INFO("angle range [%f, %f]", angle_range[0], angle_range[1]);

  Real tolx = degToRad(.5);
  Real tolfx = 1e-2;

  //set in lambda
  Real dz; // contact height error
  Vec3 point; //contact point (in world coords)
  Vec3 normal; //surface normal (in world coords)

  //lambda closure, requires C++11
  auto dotProductTangentNormal = [&] (const Real angle) -> Real
  {
    //convert contact angle to point
    Vec3 point_wheel = contactAngleToPoint(wheel_geom.radius, angle); //in wheel coords

    point = HT_wheel_to_world.applyTo(point_wheel); //convert to world coords

    //get vector tangent to wheel surface at contact point
    Vec3 tangent = -cos(angle)*HT_wheel_to_world.R.col(0) +
                    sin(angle)*HT_wheel_to_world.R.col(2);

    //get surface normal
    surfs.getDistance(point, dz, normal);

    Real dp = tangent.dot(normal);
    ROS_INFO("  angle = %f, point (wheel coords) = [%f %f %f], tangent.dot(normal) = %f", angle, point_wheel[0], point_wheel[1], point_wheel[2], dp);
    return dp;
  };

  //uncomment *one* of the following:
  //contact.angle = findRootBisection(angle_range[0],angle_range[1],dotProductTangentNormal,tolx,tolfx);
  contact.angle = findRootBrents(angle_range[0],angle_range[1],dotProductTangentNormal,tolx,tolfx); //faster

  contact.dz = dz;
  calcHTContactToWorld(point, normal, HT_wheel_to_world.R.col(1), contact.HT_world);
  contact.HT_wheel = HT_wheel_to_world.composeInvWith(contact.HT_world); //TODO, check this

  ROS_INFO("contact.angle = %f, dz = %f", contact.angle, contact.dz); //DEBUGGING
  ROS_INFO_STREAM("contact HT_wheel = \n" << contact.HT_wheel);
  ROS_INFO_STREAM("contact HT_world = \n" << contact.HT_world);
}

} //namespace
