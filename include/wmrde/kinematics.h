#ifndef _WMRDE_KINEMATICS_H_
#define _WMRDE_KINEMATICS_H_

#include <wmrde/wmrstate.h>
#include <wmrde/surface/surface.h>
#include <wmrde/contactframe.h>

namespace wmrde
{

/*!
 * calculate Jacobian of contact point velocity (i.e. velocity of point on wheel in world coords)
 * with respect to joint space velocity.
 * \param[in] mdl the WmrModel
 * \param[in] HT_world vector of transforms from each frame in model to world frame
 * \param[in] contacts vector of contact geometry for each wheel.
 * \param[out] wheel_jacobians Jacobian matrix of for all wheels in contact
 *             matrix size is 3*(num wheels in contact) x joint space size
 */
void calcWheelJacobians(
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_world,
    const std::vector<ContactFrame>& contacts,
    Matd& wheel_jacobians);

/*!
 * Compute forward velocity kinematics for a WmrModel. i.e. given actuated joint rates,
 * compute the velocity of the vehicle body frame and passive joints.
 * \param[in] mdl
 * \param[in] state
 * \param[in] u the actuated joint rates. The size must equal mdl.numActuatedFrames()
 * \param[in] contacts contact frames for all wheels in contact with ground
 * \param[in] HT_world transforms for each WmrModel frame that convert to world coordinates.
 *            This is redundant with state but is required to be precomputed. //TODO, make this optional?
 * \param[out] joint_space_vel the joint space velocity
 */
void forwardVelocityKinematics(
    const WmrModel& mdl,
    const WmrState& state,
    const Vecd& u,
    const std::vector<ContactFrame>& contacts,
    const std::vector<HTransform>& HT_world,
    JointSpaceVel& joint_space_vel);

} //namespace

#endif
