#ifndef _WMRDE_WMRSTATE_H_
#define _WMRDE_WMRSTATE_H_

#include <iostream>
#include <wmrde/wmrmodel.h>
#include <wmrde/algebra/transform.h>
#include <wmrde/algebra/rotation.h>
#include <wmrde/algebra/dynamic_matrix.h>

namespace wmrde
{

class WmrState
{
 public:
  Vec3 position; //position of body frame
  Quaternion orientation; //orientation of body frame
  Vecd joint_disp; //joint displacement vector. size is num joints (WmrModel::numJoints())

  WmrState() {}
  WmrState(
      const Vec3& body_position,
      const Quaternion body_orientation,
      const Vecd& joint_displacements)
  :
    position(body_position),
    orientation(body_orientation),
    joint_disp(joint_displacements)
  {}

  void setIdentity()
  {
    position = Vec3::Zero();
    orientation = Quaternion::Identity();
    joint_disp.setZero();
  }

  /*!
   * concatenate state into an Eigen Dynamic vector =
   * [position (3x1),
   *  orientation (4x1 quaternion coefficients),
   *  joint displacements (num joints x 1)]
   */
  Vecd toVecd() const
  {
    Vecd out;
    out.resize(7+joint_disp.rows());
    out << position, orientation.coeffs(), joint_disp;
    return out;
  }
  void fromVecd(const Vecd& in)
  {
    position = in.topRows(3); //the first 3 rows
    Eigen::Matrix<Real,4,1> tmp = in.block(3,0,4,1);
    orientation = Quaternion(tmp); //the next 4 rows
    int nj = in.rows()-7;
    joint_disp.resize(nj);
    joint_disp = in.bottomRows(nj); //the remaining rows
  }

  /*!
   * convert to state vector in which orientation is represented by
   * Euler angles instead of a quaternion. Prefer to use the quaternion
   * representation because it is less susceptible to singularities.
   */
  Vecd toVecdEuler() const
  {
    Vec3 euler;
    rotToEuler(orientation.matrix(), euler[0], euler[1], euler[2]);
    Vecd out;
    out.resize(6+joint_disp.rows());
    out << position, euler, joint_disp;
    return out;
  }
  void fromVecdEuler(const Vecd& in)
  {
    position = in.topRows(3); //the first 3 rows
    Vec3 euler = in.block(3,0,3,1); //the next 3 rows
    orientation = Quaternion(eulerToRot(euler[0], euler[1], euler[2]));
    int nj = in.rows()-6;
    joint_disp.resize(nj);
    joint_disp = in.bottomRows(nj); //the remaining rows
  }

  friend std::ostream &operator<<( std::ostream &output, const WmrState &state)
  {
    output << state.toVecdEuler();
    return output;
  }
};

/*!
 * Convert WmrState to vector of homogeneous transforms that represent the pose
 * of each frame with respect to its parent frame.
 * \param HT_parent the pose of every frame in the wmr model wrt its parent
 */
void stateToHTparent(
    const WmrModel& mdl,
    const WmrState& state,
    std::vector<HTransform>& HT_parent);

/*!
 * Given a vector of transforms for each WmrModel frame that convert to
 * parent frame coordinates, compute transforms that convert to world coordinates.
 * \param HT_parent transforms to parent frame coords
 * \param HT_world transforms to world frame coords
 */
void HTparentToHTworld(
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_parent,
    std::vector<HTransform>& HT_world);

/*!
 * 'joint space velocity' is a vector that comprises spatial velocity
 * of the body frame (in body coords) and joint rates:
 * [[wx, wy, wz, vx, vy, vz], joint_rates]
 * The total size is 6 + num joints, equal to number of degrees of freedom
 */
typedef Vecd JointSpaceVel;
typedef Vecd JointSpaceAccel; //!< the time derivative of JointSpaceVel

/*!
 * Step WmrState forward in time given joint space velocity
 * \param state the WmrState before step
 * \param joint_space_vel. 'joint space velocity'
 * \param dt the time step size
 * \return the WmrState after step
 */
WmrState stepWmrState(
    const WmrState& state,
    const JointSpaceVel& joint_space_vel,
    const Real dt);

} //namespace


#endif //_WMRDE_STATE_H_
