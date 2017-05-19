#ifndef _WMRDE_STATE_H_
#define _WMRDE_STATE_H_

#include <iostream>
#include <wmrde/wmrmodel.h>
#include <wmrde/algebra/transform.h>
#include <wmrde/algebra/rotation.h>

namespace wmrde
{

//TODO, move these typedefs?
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> Vecd;
typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> Matd;

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
  friend std::ostream &operator<<( std::ostream &output, const WmrState &state)
  {
    output << state.toVecd();
    return output;
  }
};

/*!
 * Convert WmrState to vector of homogeneous transforms which
 * represent the pose of each frame with respect to its parent frame
 * \param HT_parent the pose of every frame in the wmr model wrt its parent
 */
void stateToHTparent(
    const WmrModel& mdl,
    const WmrState& state,
    std::vector<HTransform>& HT_parent);

/*!
 * Convert vector of transforms with respect to parent frame
 * to transforms with respect to world frame
 * \param HT_parent transform to parent frame coords for every frame
 * \param HT_world transform to world frame coords
 */
void HTparentToHTworld(
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_parent,
    std::vector<HTransform>& HT_world);

/*!
 * Step WmrState forward in time given joint space velocity
 * \param state the WmrState before step
 * \param joint_space_vel. 'joint space velocity' comprises spatial velocity
 *        of the body frame (in body coords) and joint rates:
 *        [wx wy wz vx vy vz joint_rate[0] ... joint_rate[num joints-1]]
 * \param dt the time step size
 * \return the WmrState after step
 */
WmrState stepWmrState(
    const WmrState& state,
    const Vecd& joint_space_vel,
    const Real dt);

} //namespace


#endif //_WMRDE_STATE_H_
