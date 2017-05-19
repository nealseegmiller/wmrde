#ifndef _WMRDE_WMRDE_ROS_INTERFACE_H_
#define _WMRDE_WMRDE_ROS_INTERFACE_H_

#include <ros/ros.h>
#include <tf/tf.h>
#include <tf/transform_broadcaster.h>
#include <wmrde/surface/plane_surface.h>
#include <wmrde/surface/grid_surface.h>
#include <wmrde/contactgeom.h>
#include <wmrde/wmrmodel.h>
#include <visualization_msgs/Marker.h>

namespace wmrde
{

inline Mat3 TfMatrix3x3ToMat3(const tf::Matrix3x3& src)
{
  //tf::Matrix3x3 is stored internally as array of 3 Vector3's
  //one for each row

  Mat3 out;
  out << src[0].x(), src[0].y(), src[0].z(),
         src[1].x(), src[1].y(), src[1].z(),
         src[2].x(), src[2].y(), src[2].z();
  return out;
}

inline tf::Matrix3x3 Mat3ToTfMatrix3x3(const Mat3& src)
{
  return tf::Matrix3x3(
      src(0,0), src(0,1), src(0,2),
      src(1,0), src(1,1), src(1,2),
      src(2,0), src(2,1), src(2,2));
}

inline HTransform TfTransformToHTransform(const tf::Transform& src)
{
  HTransform out;
  out.R = TfMatrix3x3ToMat3(tf::Matrix3x3(src.getRotation()));
  out.t << src.getOrigin().x(), src.getOrigin().y(), src.getOrigin().z();
  return out;
}

inline tf::Transform HTransformToTfTransform(const HTransform& src)
{
  tf::Transform out;
  tf::Quaternion q;
  Mat3ToTfMatrix3x3(src.R).getRotation(q);
  out.setRotation(q);
  out.setOrigin(tf::Vector3(src.t[0], src.t[1], src.t[2]));
  return out;
}

/*!
 * broadcast the transform tree specified by HT_parent and WmrModel
 * \param tf_broadcaster non-const because .sendTransform() is non-const
 * \param mdl Specifies the names of parent and child frames for each transform
 * \param HT_parent vector of transforms from each WmrModel frame to its parent frame
 * \param stamp
 */
void tfBroadcastWmrModelTransforms(
    tf::TransformBroadcaster& tf_broadcaster,
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_parent,
    const ros::Time& stamp = ros::Time::now());

inline geometry_msgs::Point geomPoint(const double x, const double y, const double z)
{
  geometry_msgs::Point pt;
  pt.x = x;
  pt.y = y;
  pt.z = z;
  return pt;
}

inline std_msgs::ColorRGBA makeColorRGBA(const double r, const double g, const double b)
{
  std_msgs::ColorRGBA color;
  color.r = r;
  color.g = g;
  color.b = b;
  color.a = 1.0;
  return color;
}

void WheelGeomRadiusToMarker(
    const WheelGeom& wheel_geom,
    visualization_msgs::Marker& marker);

class WmrdeRosInterface // a friend class of wmrde classes so it can access private members
{
 public:
  static void PlaneSurfaceToMarker(
      const PlaneSurface& surf,
      visualization_msgs::Marker& marker);

  static void GridSurfaceToMarker(
      const GridSurface& surf,
      visualization_msgs::Marker& marker); //set marker type to either LINE_LIST or POINTS
};

} //namespace

#endif
