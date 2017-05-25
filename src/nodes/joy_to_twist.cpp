#include <ros/ros.h>
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

class JoyToTwist
{
 public:
  JoyToTwist();
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

  ros::Publisher twist_pub_;
  ros::Subscriber joy_sub_;

  double update_rate_hz_;
  ros::Timer update_timer_;
  void timerCallback(const ros::TimerEvent& ev);
};

inline double degToRad(double val) { return val*M_PI/180.0; }

JoyToTwist::JoyToTwist()
{
  pnh_ = ros::NodeHandle("~");

  //default buttons and axes based on ps3 controller
  pnh_.param("axis_linear_x", axis_linear_x_, int(Ps3::left_stick_vertical));
  pnh_.param("axis_linear_y", axis_linear_y_, int(Ps3::left_stick_horizontal));
  pnh_.param("axis_linear_z", axis_linear_z_, int(Ps3::right_stick_vertical));
  pnh_.param("axis_angular_z", axis_angular_z_, int(Ps3::right_stick_horizontal));
  pnh_.param("scale_linear", scale_linear_, 1.0);
  pnh_.param("scale_angular", scale_angular_, 1.0);

  twist_pub_ = nh_.advertise<geometry_msgs::TwistStamped>("cmd_twist", 2);
  joy_sub_ = nh_.subscribe<sensor_msgs::Joy>("joy", 2, &JoyToTwist::joyCallback, this);
}

void JoyToTwist::joyCallback(const sensor_msgs::Joy::ConstPtr& joy)
{
  geometry_msgs::Twist twist;
  twist.linear.x = scale_linear_* joy->axes[axis_linear_x_];
  twist.linear.y = scale_linear_* joy->axes[axis_linear_y_];
  twist.linear.z = scale_linear_* joy->axes[axis_linear_z_];
  twist.angular.x = 0.0;
  twist.angular.y = 0.0;
  twist.angular.z = scale_angular_* joy->axes[axis_angular_z_];

  geometry_msgs::TwistStamped twist_stamped;
  twist_stamped.header.stamp = ros::Time::now();
  twist_stamped.twist = twist;
  twist_pub_.publish(twist_stamped);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "joy_to_twist");
  JoyToTwist joy_to_twist;
  ros::spin();
}
