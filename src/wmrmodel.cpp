#include <wmrde/wmrmodel.h>
#include <wmrde/util/rosout.h>

namespace wmrde
{

Mat3 calcMomentOfInertiaBox(
    const Real mass,
    const Real lx,
    const Real ly,
    const Real lz)
{
  Mat3 out;
  out.setZero();
  Real lx2 = lx*lx;
  Real ly2 = ly*ly;
  Real lz2 = lz*lz;
  out(0,0) = mass*(ly2 + lz2)/12.0;
  out(1,1) = mass*(lx2 + lz2)/12.0;
  out(2,2) = mass*(lx2 + ly2)/12.0;
  return out;
}

Mat3 calcMomentOfInertiaCylinder(
    const Real mass,
    const Real radius,
    const Real height,
    const int axis)
{
  Real val_not_axis = mass*(3.0*radius*radius + height*height)/12.0;
  Real val_axis = mass*radius*radius/2.0;
  Mat3 out;
  for (int i = 0; i < 3; i++)
  {
    out(i,i) = (i == axis) ? val_axis : val_not_axis;
  }
}

int WmrModel::frameNameToIndex(const std::string& name) const
{
  for (int i = 0; i < numFrames(); i++)
  {
    if (name.compare(frames_[i].name) == 0) { return i; }
  }
  return -1;
}

bool WmrModel::addBodyFrame(const std::string& name)
{
  if (!frames_.empty())
  {
    ROS_ERROR("Failed to add body frame, frames already exist.");
    return false;
  }
  if (name.empty())
  {
    ROS_ERROR("Failed to add body frame, name is empty()");
    return false;
  }

  Frame frame(name, Frame::FREE, -1, HTransform::Identity(), false, false);
  frames_.push_back(std::move(frame));
  return true;
}

bool WmrModel::addFrame(const Frame& frame)
{
  if (frames_.empty())
  {
    ROS_ERROR("Failed to add frame, body frame must be added first.");
    return false;
  }
  //check if frame is valid
  if (frame.name.empty())
  {
    ROS_ERROR("Failed to add frame, name is empty()");
    return false;
  }
  if (frameNameToIndex(frame.name) > 0)
  {
    ROS_ERROR("Failed to add frame, name '%s' is already in use.", frame.name.c_str());
    return false;
  }
  if (frame.dof_type < 0 || frame.dof_type > 5)
  {
    ROS_ERROR("Failed to add frame, dof_type = %d is invalid", frame.dof_type);
    return false;
  }
  if (frame.wheel_geom && frame.dof_type != Frame::REV_Y)
  {
    ROS_ERROR("Failed to add frame, wheel frame dof_type %d != %d", frame.dof_type, Frame::REV_Y);
    return false;
  }
  if (frame.parent_idx < 0 || frame.parent_idx >= numFrames())
  {
    ROS_ERROR("Failed to add frame, parent_idx = %d is invalid", frame.parent_idx);
    return false;
  }

  frames_.push_back(frame.deepcopy());

  //update member vars
  int frame_idx = numFrames()-1;
  if (frame.wheel_geom)
  {
    num_wheel_frames_++;
    wheel_frame_inds_.push_back(frame_idx);
  }
  if (frame.is_actuated)
  {
    num_actuated_frames_++;
    actuated_frame_inds_.push_back(frame_idx);
  }

  return true;
}

} //namespace
