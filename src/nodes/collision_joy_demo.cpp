#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>

#include <wmrde/algebra/rotation.h>
#include <wmrde/surface/plane_surface.h>
#include <wmrde/surface/grid_surface.h>
#include <wmrde/wmrde_ros_interface.h>
#include <wmrde/collision.h>

#include <visualization_msgs/MarkerArray.h>

namespace vm = visualization_msgs;

using namespace wmrde;

class CollisionJoyDemo
{
 public:
    CollisionJoyDemo();

 private:

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;

  ros::Publisher surfs_marker_pub_;
  ros::Publisher wheel_marker_pub_;
  ros::Publisher distance_pub_;

  //params
  std::string world_frame_;
  std::string wheel_frame_;
  std::string contact_frame_;
  bool wheel_contact_frames_coincident_; //if true set the contact frame to be coincident with the wheel frame,
  //this allows testing of Surfaces class in isolation
  Surfaces surfs_;
  WheelGeom wheel_geom_;

  tf::TransformListener tf_listener_;
  tf::TransformBroadcaster tf_broadcaster_;

  double update_rate_hz_;
  ros::Timer update_timer_;
  void timerCallback(const ros::TimerEvent& ev);

  void initSurfaces();
};

void CollisionJoyDemo::initSurfaces()
{
  surfs_.clear();

  //get parameters
  bool enable_plane, enable_grid;
  std::string grid_surface_filename;

  pnh_.param("enable_plane", enable_plane, true);
  pnh_.param("enable_grid", enable_grid, true);
  pnh_.param("grid_surface_filename", grid_surface_filename, std::string());

  //init the surfaces
  PlaneSurface plane_surf;
  if (enable_plane)
  {
    HTransform T;
    T.setIdentity();
    T.t << 0.0, 0.0, 0.25;
    T.R = eulerToRot(degToRad(0.0), degToRad(10.0), degToRad(0.0));
    BoundingBox bb(0.0, 1.0, 0.0, 2.0);
    plane_surf = PlaneSurface(T, bb);
    surfs_.addSurface(plane_surf);
  }

  GridSurface grid_surf;
  if (enable_grid)
  {
    grid_surf = GridSurface(grid_surface_filename);
    surfs_.addSurface(grid_surf);
  }

  //publish marker visualization for surfaces
  {
    vm::MarkerArray ma;

    vm::Marker m_init;
    m_init.header.frame_id = world_frame_;
    m_init.header.stamp = ros::Time::now();

    if (plane_surf.initialized())
    {
      vm::Marker marker = m_init;
      marker.id = 0;
      marker.color = makeColorRGBA(0.0, 0.0, 1.0);
      marker.scale.x = 0.02;
      WmrdeRosInterface::PlaneSurfaceToMarker(plane_surf, marker);
      ma.markers.push_back(marker);
    }

    if (grid_surf.initialized())
    {
      vm::Marker marker = m_init;
      marker.type = vm::Marker::POINTS;
      marker.id = 1;
      marker.color = makeColorRGBA(0.0, 1.0, 0.0);
      marker.scale.x = 0.01;
      WmrdeRosInterface::GridSurfaceToMarker(grid_surf, marker);
      ma.markers.push_back(marker);
    }

    surfs_marker_pub_.publish(ma);
  }
}

CollisionJoyDemo::CollisionJoyDemo()
{
  ros::Duration(1.0).sleep(); //sleep so no rosout is missed
  pnh_ = ros::NodeHandle("~");

  //get parameters
  pnh_.param("world_frame", world_frame_, std::string("world"));
  pnh_.param("wheel_frame", wheel_frame_, std::string("wheel"));
  pnh_.param("contact_frame", contact_frame_, std::string("contact"));
  pnh_.param("wheel_contact_frames_coincident", wheel_contact_frames_coincident_, false);
  pnh_.param("wheel_radius", wheel_geom_.radius, 0.3);

  //init publisher
  surfs_marker_pub_ = pnh_.advertise<vm::MarkerArray>("surface_marker", 2, true);
  wheel_marker_pub_ = pnh_.advertise<vm::Marker>("wheel_marker", 2);
  distance_pub_ = pnh_.advertise<vm::MarkerArray>("distance_to_surface", 2);

  initSurfaces();

  pnh_.param("update_rate_hz", update_rate_hz_, 10.0);
  update_timer_ = nh_.createTimer(ros::Duration(1.0/update_rate_hz_), &CollisionJoyDemo::timerCallback, this);

}

void CollisionJoyDemo::timerCallback(const ros::TimerEvent& ev)
{
  //get point from tf listener
  tf::StampedTransform tf_wheel_to_world;
  try
  {
    tf_listener_.lookupTransform(world_frame_, wheel_frame_, ros::Time(0), tf_wheel_to_world);
  }
  catch (tf::TransformException& ex)
  {
    ROS_ERROR("%s", ex.what());
    return;
  }
  HTransform HT_wheel_to_world = TfTransformToHTransform(tf::Transform(tf_wheel_to_world));

  //publish marker for wheel frame
  {
    vm::Marker marker;
    WheelGeomRadiusToMarker(wheel_geom_, marker);
    marker.header.frame_id = wheel_frame_;
    marker.header.stamp = ros::Time::now();
    marker.id = 0;
    marker.color = makeColorRGBA(0.5, 0.5, 0.5);
    marker.scale.x = 0.01; //only x is used

    wheel_marker_pub_.publish(marker);
  }

  //get the wheel-ground contact frame
  ContactGeom contact;
  if (wheel_contact_frames_coincident_)
  {
    contact.HT_wheel.setIdentity();
    contact.HT_world = HT_wheel_to_world;
  }
  else
  {
    calcContactGeomRootfinding(surfs_, wheel_geom_, HT_wheel_to_world, contact);
  }

  //broadcast the contact frame
  tf::Transform tf_contact_to_wheel = HTransformToTfTransform(contact.HT_wheel);
  tf_broadcaster_.sendTransform(tf::StampedTransform(tf_contact_to_wheel, ros::Time::now(), wheel_frame_, contact_frame_));

  //get distance of contact point from surface
  Vec3 point  = contact.HT_world.t;
  ROS_DEBUG("contact point = %f, %f, %f", point[0], point[1], point[2]);

  //init marker
  vm::MarkerArray ma;
  vm::Marker m_init;
  m_init.header.frame_id = world_frame_;
  m_init.header.stamp = ros::Time::now();
  m_init.type = vm::Marker::ARROW;
  m_init.scale.x = 0.01;
  m_init.scale.y = 0.02;
  m_init.scale.z = 0.0;

  //pub arrow marker from point to nearest point on surface
  {
    Vec3 normal;
    Real distance;
    surfs_.getDistance(point, distance, normal);

    Vec3 nearest = point - distance*normal;
    ROS_DEBUG("distance = %f, normal = %f, %f, %f", distance, normal[0], normal[1], normal[2]);
    ROS_DEBUG("nearest point = %f, %f, %f", nearest[0], nearest[1], nearest[2]);

    vm::Marker marker = m_init;
    marker.id = 0;
    marker.color = makeColorRGBA(0.0, 1.0, 0.0);
    marker.points = {
        geomPoint(point[0], point[1], point[2]),
        geomPoint(nearest[0], nearest[1], nearest[2])
    };
    ma.markers.push_back(marker);
  }

  //pub another arrow from point to point on point directly below on surface
  {
    double height;
    surfs_.getHeight(point, height);
    Vec3 below = point;
    below[2] = height;

    ROS_DEBUG("height = %f", height);
    ROS_DEBUG("below pt = %f, %f, %f\n", below[0], below[1], below[2]);

    vm::Marker marker = m_init;
    marker.id = 1;
    marker.color = makeColorRGBA(0.0, 1.0, 1.0);
    marker.points = {
        geomPoint(point[0], point[1], point[2]),
        geomPoint(below[0], below[1], below[2])
    };
    ma.markers.push_back(marker);
  }

  distance_pub_.publish(ma);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "surface_joy_demo");

  CollisionJoyDemo demo;
  ros::spin();
}
