#include <tf/transform_broadcaster.h>

#include <wmrde/algebra/rotation.h>
#include <wmrde/surface/plane_surface.h>
#include <wmrde/surface/grid_surface.h>
#include <wmrde/wmrde_ros_interface.h>
#include <wmrde/collision.h>

#include <geometry_msgs/TwistStamped.h>
#include <visualization_msgs/MarkerArray.h>

namespace vm = visualization_msgs;

using namespace wmrde;

class CollisionDemo
{
 public:
  CollisionDemo();

 private:

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;

  ros::Subscriber twist_sub_;
  ros::Publisher surfs_marker_pub_;
  ros::Publisher wheel_marker_pub_;
  ros::Publisher distance_pub_;

  //params
  bool teleop_contact_frame_directly_;

  geometry_msgs::Twist twist_;
  HTransform HT_wheel_to_world_; //transform from wheel frame to world
  Surfaces surfs_;
  WheelGeom wheel_geom_;

  tf::TransformBroadcaster tf_broadcaster_;

  double update_rate_hz_;
  ros::Timer update_timer_;
  void timerCallback(const ros::TimerEvent& ev);

  void initTransform();
  void initSurfaces();
  void twistCallback(const geometry_msgs::TwistStampedConstPtr& msg); //TODO
};

void CollisionDemo::initTransform()
{
  Real init_x, init_y, init_z, init_roll_deg, init_pitch_deg, init_yaw_deg;
  pnh_.param("init_x", init_x, 0.0);
  pnh_.param("init_y", init_y, 0.0);
  pnh_.param("init_z", init_z, 0.5);
  pnh_.param("init_roll_deg", init_roll_deg, 0.0);
  pnh_.param("init_pitch_deg", init_pitch_deg, 0.0);
  pnh_.param("init_yaw_deg", init_yaw_deg, 0.0);

  HT_wheel_to_world_.setIdentity();
  HT_wheel_to_world_.t << init_x, init_y, init_z;
  HT_wheel_to_world_.R = eulerToRot(
      degToRad(init_roll_deg),
      degToRad(init_pitch_deg),
      degToRad(init_yaw_deg));

}

void CollisionDemo::initSurfaces()
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
    BoundingBox bb(-0.5, 0.5, -1.0, 1.0);
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
    m_init.header.frame_id = "world";
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

CollisionDemo::CollisionDemo()
{
  ros::Duration(1.0).sleep(); //sleep so no rosout is missed
  pnh_ = ros::NodeHandle("~");

  //get parameters
  //Ordinarily user teleops the wheel frame with twist commands, and contact frame
  //is calculated by collision check with ground surfaces. but option is provided to
  //teleop the contact frame directly so surface classes can be tested in isolation.
  pnh_.param("teleop_contact_frame_directly", teleop_contact_frame_directly_, false);
  pnh_.param("wheel_radius", wheel_geom_.radius, 0.3);

  //init publisher
  twist_sub_ = nh_.subscribe<geometry_msgs::TwistStamped>("cmd_twist", 2, &CollisionDemo::twistCallback, this);
  surfs_marker_pub_ = pnh_.advertise<vm::MarkerArray>("surface_marker", 2, true);
  wheel_marker_pub_ = pnh_.advertise<vm::Marker>("wheel_marker", 2);
  distance_pub_ = pnh_.advertise<vm::MarkerArray>("distance_to_surface", 2);

  initTransform();
  initSurfaces();

  pnh_.param("update_rate_hz", update_rate_hz_, 10.0);
  update_timer_ = nh_.createTimer(ros::Duration(1.0/update_rate_hz_), &CollisionDemo::timerCallback, this);

  void twistCallback(const geometry_msgs::TwistStampedConstPtr& msg); //TODO
}

void CollisionDemo::twistCallback(const geometry_msgs::TwistStampedConstPtr& msg)
{
  twist_ = msg->twist;
}

void CollisionDemo::timerCallback(const ros::TimerEvent& ev)
{
  //update teleop frame according to twist
  Real dt = 1.0/update_rate_hz_;

  //update translation
  Vec3 linear_vel;
  linear_vel << twist_.linear.x, twist_.linear.y, twist_.linear.z;
  HT_wheel_to_world_.t += HT_wheel_to_world_.R*linear_vel*dt;

  //update rotation
  Quaternion q(HT_wheel_to_world_.R);
  Vec3 angular_vel;
  angular_vel << twist_.angular.x, twist_.angular.y, twist_.angular.z;
  q = stepRotationQuat(q, angular_vel, dt);
  HT_wheel_to_world_.R = q.matrix();

//  ROS_INFO_STREAM("HT_wheel_to_world = \n" << HT_wheel_to_world_);

  Vec3 point; //contact point
  if (!teleop_contact_frame_directly_)
  {
    ContactFrame contact;
    calcContactFrameRootfind(surfs_, wheel_geom_, HT_wheel_to_world_, contact);

    //publish marker for wheel frame
    {
      vm::Marker marker;
      WheelGeomRadiusToMarker(wheel_geom_, marker);
      marker.header.frame_id = "wheel";
      marker.header.stamp = ros::Time::now();
      marker.id = 0;
      marker.color = makeColorRGBA(0.5, 0.5, 0.5);
      marker.scale.x = 0.01; //only x is used

      wheel_marker_pub_.publish(marker);
    }

    //broadcast transforms
    tf::Transform transform = HTransformToTfTransform(HT_wheel_to_world_);
    tf_broadcaster_.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "wheel"));

    transform = HTransformToTfTransform(contact.HT_wheel);
    tf_broadcaster_.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "wheel", "contact"));

    point = contact.HT_world.t;
  }
  else //wheel frame is contact frame
  {
    //broadcast transforms
    tf::Transform transform = HTransformToTfTransform(HT_wheel_to_world_);
    tf_broadcaster_.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "contact"));

    point = HT_wheel_to_world_.t;
  }

  ROS_INFO_STREAM("contact point = " << point.transpose());

  //init marker
  vm::MarkerArray ma;
  vm::Marker m_init;
  m_init.header.frame_id = "world";
  m_init.header.stamp = ros::Time::now();
  m_init.type = vm::Marker::ARROW;
  m_init.scale.x = 0.01;
  m_init.scale.y = 0.02;
  m_init.scale.z = 0.0;

  //pub arrow marker from point to nearest point on surface
  {
    Vec3 normal;
    Real dz; //contact height error
    surfs_.getDistance(point, dz, normal);

    Vec3 nearest = point - dz*normal;
    ROS_INFO_STREAM("dz = " << dz << ", normal = " << normal.transpose());
    ROS_INFO_STREAM("nearest point = " << nearest.transpose());

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

    ROS_INFO("surface height = %f\n", height);

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
  ros::init(argc, argv, "collision_demo");

  CollisionDemo demo;
  ros::spin();
}
