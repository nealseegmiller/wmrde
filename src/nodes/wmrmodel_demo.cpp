#include <ros/ros.h>
#include <tf/transform_broadcaster.h>

#include <wmrde/demo/zoemodel.h>
#include <wmrde/wmrde_ros_interface.h>
#include <wmrde/wmrstate.h>

#include <sensor_msgs/Joy.h>
#include <visualization_msgs/MarkerArray.h>


namespace vm = visualization_msgs;

using namespace wmrde;

class WmrModelDemo
{
 public:
  WmrModelDemo();

 private:

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;

  ros::Subscriber joy_sub_;
  ros::Publisher mesh_pub_;

  //params
  int joy_button_idx_; //index of joy button to toggle through joints
  int joy_axis_idx_; //index of joy axis to actuate joint
  double scale_joint_vel_; //multiply axis value by this scale to get velocity

  WmrModel mdl_;
  vm::MarkerArray mdl_markers_;
  WmrState state_;
  int joint_idx_; //active joint idx, (i.e. element of joint space velocity)
  double joint_vel_;
  bool button_pressed_;

  tf::TransformBroadcaster tf_broadcaster_;

  double update_rate_hz_;
  ros::Timer update_timer_;
  void timerCallback(const ros::TimerEvent& ev);

  void joyCallback(const sensor_msgs::JoyConstPtr& joy);
};


WmrModelDemo::WmrModelDemo()
:
  joint_idx_(0),
  joint_vel_(0.0),
  button_pressed_(false)
{
  ros::Duration(1.0).sleep(); //sleep so no rosout is missed
  pnh_ = ros::NodeHandle("~");

  //get parameters
  pnh_.param("joy_button_idx", joy_button_idx_, 12); //ps3 triangle button
  pnh_.param("joy_axis_idx", joy_axis_idx_, 3); //ps3 right stick vertical
  pnh_.param("scale_joint_vel", scale_joint_vel_, 1.0);

  //init publisher/subscriber
  joy_sub_ = nh_.subscribe<sensor_msgs::Joy>("joy", 2, &WmrModelDemo::joyCallback, this);
  mesh_pub_ = pnh_.advertise<vm::MarkerArray>("wmrmodel_mesh", 2, true); //latched

  //init the model and marker array
  makeZoeModel(mdl_);
  makeZoeModelMarkers(mdl_, mdl_markers_);

  //init the state
  state_.joint_disp.resize(mdl_.numJoints());
  state_.setIdentity();

  ROS_INFO_STREAM("initial state = \n" << state_);

  pnh_.param("update_rate_hz", update_rate_hz_, 10.0);
  update_timer_ = nh_.createTimer(ros::Duration(1.0/update_rate_hz_), &WmrModelDemo::timerCallback, this);
}

void WmrModelDemo::joyCallback(const sensor_msgs::JoyConstPtr& joy)
{
  if (!button_pressed_ && joy->buttons[joy_button_idx_])
  {
    //button changed from not pressed to pressed
    joint_idx_++;
    joint_idx_ = joint_idx_ % (mdl_.numDof()); //wrap around back to zero if exceed size of joint space
  }
  button_pressed_ = joy->buttons[joy_button_idx_];

  joint_vel_ = scale_joint_vel_*joy->axes[joy_axis_idx_];
}

void WmrModelDemo::timerCallback(const ros::TimerEvent& ev)
{
  //TODO, set joint space velocity according to joint_idx_, joint_vel_ and update state
  ROS_INFO("joint_idx = %d, joint_vel = %f", joint_idx_, joint_vel_);
  if (std::abs(joint_vel_) > 0)
  {
    Vecd qvel;
    qvel.resize(mdl_.numDof());
    qvel.setZero();
    qvel[joint_idx_] = joint_vel_;
    state_ = stepWmrState(state_, qvel, 1.0/update_rate_hz_);
    ROS_INFO_STREAM("qvel = \n" << qvel);
    ROS_INFO_STREAM("state (after step) = \n" << state_);
  }

  //compute HT_parent
  std::vector<HTransform> HT_parent;
  stateToHTparent(mdl_, state_, HT_parent);

  //broadcast transforms over tf
  tfBroadcastWmrModelTransforms(tf_broadcaster_, mdl_, HT_parent, ros::Time::now());

  //TODO, must republish marker array to update render in rviz?
  if (!mdl_markers_.markers.empty())
  {
    vm::MarkerArray ma = mdl_markers_;
    //set stamps
    ros::Time cur_time = ros::Time::now();
    for (vm::Marker& marker : ma.markers) { marker.header.stamp = cur_time; }
    mesh_pub_.publish(ma);
    ROS_INFO("published MarkerArray with %zu MESH_RESOURCE markers\n", ma.markers.size());
  }

}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "wmrmodel_demo");

  WmrModelDemo demo;
  ros::spin();
}

