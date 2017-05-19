#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/Joy.h>

//button and axis indices for Ps3 controller
struct Ps3
{
  enum button_indices {
    dpad_up_button = 4,
    dpad_right_button = 5,
    dpad_down_button = 6,
    dpad_left_button = 7,
    triangle_button = 12,
    circle_button = 13,
    x_button = 14,
    square_button = 15,
  };
  enum axis_indices {
    left_stick_horizontal = 0,
    left_stick_vertical = 1,
    right_stick_horizontal = 2,
    right_stick_vertical = 3,
  };
};

class JoyTeleop
{
 public:
  JoyTeleop();
  void initTransform();
  double rescaleAxis(double val) const ;

 private:

  void joyCallback(const sensor_msgs::Joy::ConstPtr& joy);

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;

  int axis_linear_x_;
  int axis_linear_y_;
  int axis_linear_z_;
  int axis_angular_z_;
  double scale_linear_;
  double scale_angular_;
  int init_transform_button_;

  ros::Publisher twist_pub_;
  ros::Subscriber joy_sub_;
  geometry_msgs::Twist twist_;

  bool enable_tf_; //broadcast a tf that moves according to cmd_twist

  std::string parent_frame_;
  std::string child_frame_;
  tf::TransformBroadcaster tf_broadcaster_;
  tf::Transform transform_;

  double update_rate_hz_;
  ros::Timer update_timer_;
  void timerCallback(const ros::TimerEvent& ev);
};

inline double degToRad(double val) { return val*M_PI/180.0; }

void JoyTeleop::initTransform()
{
  ROS_INFO("Initialize transform.");
  double init_x, init_y, init_z;
  pnh_.param("init_x", init_x, 0.0);
  pnh_.param("init_y", init_y, 0.0);
  pnh_.param("init_z", init_z, 0.0);
  transform_.setOrigin(tf::Vector3(init_x, init_y, init_z));

  double init_roll, init_pitch, init_yaw;
  pnh_.param("init_roll_deg", init_roll, 0.0);
  pnh_.param("init_pitch_deg", init_pitch, 0.0);
  pnh_.param("init_yaw_deg", init_yaw, 0.0);
  transform_.setRotation(tf::createQuaternionFromRPY(
      degToRad(init_roll),
      degToRad(init_pitch),
      degToRad(init_yaw)));
}

JoyTeleop::JoyTeleop()
{
  pnh_ = ros::NodeHandle("~");

  //default buttons and axes based on ps3 controller
  pnh_.param("axis_linear_x", axis_linear_x_, int(Ps3::left_stick_vertical));
  pnh_.param("axis_linear_y", axis_linear_y_, int(Ps3::left_stick_horizontal));
  pnh_.param("axis_linear_z", axis_linear_z_, int(Ps3::right_stick_vertical));
  pnh_.param("axis_angular_z", axis_angular_z_, int(Ps3::right_stick_horizontal));
  pnh_.param("scale_linear", scale_linear_, 1.0);
  pnh_.param("scale_angular", scale_angular_, 1.0);

  pnh_.param("init_transform_button", init_transform_button_, int(Ps3::square_button));

  twist_pub_ = nh_.advertise<geometry_msgs::TwistStamped>("cmd_twist", 2);
  joy_sub_ = nh_.subscribe<sensor_msgs::Joy>("joy", 2, &JoyTeleop::joyCallback, this);

  pnh_.param("enable_tf", enable_tf_, true);
  if (enable_tf_)
  {
    pnh_.param("parent_frame", parent_frame_, std::string("world"));
    pnh_.param("child_frame", child_frame_, std::string("teleop_frame"));
    initTransform();

    pnh_.param("update_rate_hz", update_rate_hz_, 10.0);
    update_timer_ = nh_.createTimer(ros::Duration(1.0/update_rate_hz_), &JoyTeleop::timerCallback, this);
  }
}

template <typename T>
double sign(T val)
{
  return (val > 0) - (val < 0);
}

void JoyTeleop::joyCallback(const sensor_msgs::Joy::ConstPtr& joy)
{
  twist_.linear.x = scale_linear_* joy->axes[axis_linear_x_];
  twist_.linear.y = scale_linear_* joy->axes[axis_linear_y_];
  twist_.linear.z = scale_linear_* joy->axes[axis_linear_z_];
  twist_.angular.x = 0.0;
  twist_.angular.y = 0.0;
  twist_.angular.z = scale_angular_* joy->axes[axis_angular_z_];

  geometry_msgs::TwistStamped twist_stamped;
  twist_stamped.header.stamp = ros::Time::now();
  twist_stamped.twist = twist_;
  twist_pub_.publish(twist_stamped);

  if (enable_tf_ && joy->buttons[init_transform_button_]) { initTransform(); }
}

void JoyTeleop::timerCallback(const ros::TimerEvent& ev)
{
  //update transform according to twist
  double dt = 1.0/update_rate_hz_; //time step

  //update origin
  tf::Vector3 linear(twist_.linear.x, twist_.linear.y, twist_.linear.z);
  tf::Matrix3x3 R(transform_.getRotation()); //rotation matrix
  linear = R*linear; //rotate into parent coords
  transform_.setOrigin(transform_.getOrigin() + linear*dt);

  //update rotation
  tf::Vector3 angular_vel(twist_.angular.x, twist_.angular.y, twist_.angular.z);
  double angle = angular_vel.length()*dt;
  tf::Vector3 axis = angular_vel/angular_vel.length(); //unit vector

  if (angle > 0.0)
  {
    tf::Quaternion dq;
    dq.setRotation(R*axis, angle); //differential rotation
    tf::Quaternion q = dq * transform_.getRotation(); //compose rotations
    q.normalize();
    transform_.setRotation(q);
  }
  tf_broadcaster_.sendTransform(tf::StampedTransform(transform_, ros::Time::now(), parent_frame_, child_frame_));
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "joy_teleop");
  JoyTeleop joy_teleop;
  ros::spin();
}
