#include <wmrde/wmrde_ros_interface.h>

namespace vm = visualization_msgs;

namespace wmrde
{

void tfBroadcastWmrModelTransforms(
    tf::TransformBroadcaster& tf_broadcaster, //can't be const
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_parent,
    const ros::Time& stamp,
    const std::string& world_frame)
{
  for (size_t i = 0; i < HT_parent.size(); i++)
  {
    tf::Transform tf_parent = HTransformToTfTransform(HT_parent[i]);
    int parent_idx = mdl.getFrame(i).parent_idx;
    std::string parent_frame = parent_idx >= 0 ? mdl.getFrame(parent_idx).name : world_frame;
    std::string child_frame = mdl.getFrame(i).name;
    tf_broadcaster.sendTransform(tf::StampedTransform(tf_parent, stamp, parent_frame, child_frame));
  }
}

void tfBroadcastContactTransforms(
    tf::TransformBroadcaster& tf_broadcaster,
    const WmrModel& mdl,
    const std::vector<ContactFrame>& contacts,
    const ros::Time& stamp)
{
  for (const ContactFrame& contact : contacts)
  {
    tf::Transform tf_parent = HTransformToTfTransform(contact.HT_wheel);
    int parent_idx = contact.parent_idx;
    std::string parent_frame = mdl.getFrame(contact.parent_idx).name;
    std::string child_frame = parent_frame + "_contact";
    tf_broadcaster.sendTransform(tf::StampedTransform(tf_parent, stamp, parent_frame, child_frame));
  }
}

void WheelGeomRadiusToMarker(
    const WheelGeom& wheel_geom,
    vm::Marker& marker)
{
  marker.type = vm::Marker::LINE_STRIP;
  marker.points.clear();
  int num_points = 73;
  double angle_step = 2.0*M_PI/(num_points-1);
  double r = wheel_geom.radius;
  for (int i = 0; i < num_points; i++)
  {
    double angle = i*angle_step;
    marker.points.push_back(geomPoint(cos(angle)*r, 0.0, sin(angle)*r));
  }
}

void WmrdeRosInterface::PlaneSurfaceToMarker(
    const PlaneSurface& surf,
    vm::Marker& marker)
{
  marker.type = vm::Marker::LINE_STRIP;
  marker.points.clear();

  //get bounding box
  BoundingBox bb;
  if (surf.bounded_) //can access private members because friend class
  {
    bb = surf.bb_;
  }
  else
  {
    //unit square to visualize unbounded plane
    bb.x_min = 0.0;
    bb.x_max = 1.0;
    bb.y_min = 0.0;
    bb.y_max = 1.0;
  }

  //get bb points
  std::array<Vec3,5> pts;
  pts[0] << bb.x_min, bb.y_min, 0.0;
  pts[1] << bb.x_max, bb.y_min, 0.0;
  pts[2] << bb.x_max, bb.y_max, 0.0;
  pts[3] << bb.x_min, bb.y_max, 0.0;
  pts[4] = pts[0]; //close

  //transform the points and copy to marker
  for (Vec3& pt : pts)
  {
    pt = surf.T_.applyTo(pt);
    marker.points.push_back(geomPoint(pt[0], pt[1], pt[2]));
  }
}

void WmrdeRosInterface::GridSurfaceToMarker(
    const GridSurface& surf,
    vm::Marker& marker)
{
  marker.points.clear();

  //abbreviate
  const GridVector& X = surf.G_.getX(0);
  const GridVector& Y = surf.G_.getX(1);

  //lambda to get geometry_msgs::Point from grid surface
  auto getPoint = [&surf] (const int i, const int j) -> geometry_msgs::Point
  {
    Vec3 pt = surf.getPoint(i,j);
    return geomPoint(pt[0], pt[1], pt[2]);
  };

  for (int i = 0; i < X.numel; i++)
  {
    for (int j = 0; j < Y.numel; j++)
    {
      if (marker.type == vm::Marker::LINE_LIST)
      {
        if (i+1 < X.numel)
        {
          marker.points.push_back(getPoint(i,j));
          marker.points.push_back(getPoint(i+1,j));
        }
        if (j+1 < Y.numel)
        {
          marker.points.push_back(getPoint(i,j));
          marker.points.push_back(getPoint(i,j+1));
        }
      }
      else if (marker.type == vm::Marker::POINTS)
      {
        marker.points.push_back(getPoint(i,j));
      }
      else //default to points
      {
        marker.type = vm::Marker::POINTS;
        marker.points.push_back(getPoint(i,j));
      }
    }
  }
}

} //namespace
