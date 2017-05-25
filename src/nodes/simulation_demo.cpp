#include <ros/ros.h>
#include <tf/transform_broadcaster.h>

#include <wmrde/demo/zoemodel.h>
#include <wmrde/wmrde_ros_interface.h>
#include <wmrde/init_contact.h>

#include <geometry_msgs/TwistStamped.h>
#include <visualization_msgs/MarkerArray.h>

namespace vm = visualization_msgs;

using namespace wmrde;

class SimulationDemo
{
 public:
  SimulationDemo();

 private:

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;

  ros::Subscriber twist_sub_; //TODO
  ros::Publisher wmrmodel_mesh_pub_; //for MESH_RESOURCE markers for WmrModel frames
  ros::Publisher surfs_marker_pub_;

  //params
  std::string world_frame_;

  geometry_msgs::Twist twist_;
  Surfaces surfs_;
  WmrModel mdl_;
  vm::MarkerArray mdl_markers_;
  WmrState state_;

  tf::TransformBroadcaster tf_broadcaster_;

  double update_rate_hz_;
  ros::Timer update_timer_;
  void timerCallback(const ros::TimerEvent& ev);

  void initState();
  void initSurfaces();
  void twistCallback(const geometry_msgs::TwistStampedConstPtr& msg); //TODO
};

void SimulationDemo::initState()
{
  state_.joint_disp.resize(mdl_.numJoints());
  state_.setIdentity();

  Vecd state_vec = state_.toVecdEuler();
  std::vector<double> data;
  pnh_.getParam("state_vec", data);
  for (int i = 0; i < state_vec.rows(); i++)
  {
    if (i >= (int) data.size()) { break; }
    state_vec[i] = data[i];
  }
  ROS_INFO_STREAM("initial state = \n" << state_vec);
  state_.fromVecdEuler(state_vec);
}

void SimulationDemo::initSurfaces()
{
  surfs_.clear();

  //init the surfaces
  PlaneSurface plane_surf[2];
  {
    HTransform T;
    T.setIdentity();
    T.t << 0.0, 0.0, 0.25;
    T.R = eulerToRot(degToRad(0.0), degToRad(10.0), degToRad(0.0));
    BoundingBox bb(0.0, 1.0, 0.0, 2.0);
    plane_surf[0] = PlaneSurface(T, bb);
    surfs_.addSurface(plane_surf[0]);

    T.setIdentity();
    plane_surf[1] = PlaneSurface(T);
    surfs_.addSurface(plane_surf[1]);
  }

  //publish marker array to visualize surfaces
  {
    vm::MarkerArray ma;

    vm::Marker m_init;
    m_init.header.frame_id = "world";
    m_init.header.stamp = ros::Time::now();
    m_init.color = makeColorRGBA(0.0, 0.0, 1.0);
    m_init.scale.x = 0.02;

    for (size_t i = 0; i < 2; i++)
    {
      vm::Marker marker = m_init;
      marker.id = i;
      WmrdeRosInterface::PlaneSurfaceToMarker(plane_surf[i], marker);
      ma.markers.push_back(marker);
    }

    surfs_marker_pub_.publish(ma);
  }
}

SimulationDemo::SimulationDemo()
{
  ros::Duration(1.0).sleep(); //sleep so no rosout is missed
  pnh_ = ros::NodeHandle("~");

  //get parameters

  //init publisher/subscriber
  twist_sub_ = nh_.subscribe<geometry_msgs::TwistStamped>("cmd_twist", 2, &SimulationDemo::twistCallback, this);
  wmrmodel_mesh_pub_ = pnh_.advertise<vm::MarkerArray>("wmrmodel_mesh", 2, true); //latched
  surfs_marker_pub_ = pnh_.advertise<vm::MarkerArray>("surface_marker", 2, true); //latched

  //init the model and marker array
  makeZoeModel(mdl_);
  makeZoeModelMarkers(mdl_, mdl_markers_);

  initState();
  initSurfaces();

  pnh_.param("update_rate_hz", update_rate_hz_, 10.0);
  update_timer_ = nh_.createTimer(ros::Duration(1.0/update_rate_hz_), &SimulationDemo::timerCallback, this);
}

void SimulationDemo::twistCallback(const geometry_msgs::TwistStampedConstPtr& msg)
{
  twist_ = msg->twist;
}

void SimulationDemo::timerCallback(const ros::TimerEvent& ev)
{
  //compute HT_world
  std::vector<HTransform> HT_parent, HT_world;
  stateToHTparent(mdl_, state_, HT_parent);
  HTparentToHTworld(mdl_, HT_parent, HT_world);

  //init contact with surface
  std::vector<ContactFrame> contacts;
//  calcWmrModelContactFrames(mdl_, surfs_, HT_world, contacts); //DEBUGGING
  initWmrModelSurfaceContact(mdl_, surfs_, state_, contacts);

  //broadcast transforms over tf
  ros::Time cur_time = ros::Time::now();
  tfBroadcastWmrModelTransforms(tf_broadcaster_, mdl_, HT_parent, cur_time);
  tfBroadcastContactTransforms(tf_broadcaster_, mdl_, contacts, cur_time);

  //TODO, must republish marker array to update render in rviz?
  if (!mdl_markers_.markers.empty())
  {
    vm::MarkerArray ma = mdl_markers_;
    //set stamps
    ros::Time cur_time = ros::Time::now();
    for (vm::Marker& marker : ma.markers) { marker.header.stamp = cur_time; }
    wmrmodel_mesh_pub_.publish(ma);
    ROS_INFO("published MarkerArray with %zu MESH_RESOURCE markers\n", ma.markers.size());
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "wmrmodel_demo");

  SimulationDemo demo;
  ros::spin();
}

