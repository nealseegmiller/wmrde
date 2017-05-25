#include <wmrde/init_contact.h>
#include <wmrde/kinematics.h>
#include <wmrde/util/index_util.h>
#include <wmrde/util/linesearch.h>
#include <wmrde/util/rosout.h>

namespace wmrde
{

//TODO, move this?
/*!
 * convert joint space velocity into time derivative of state (in which
 * orientation is expressed using Euler angles)
 * \param R_world rotation matrix from body to world coords
 * \param roll
 * \param pitch
 * \param V each column is a joint space velocity vector
 * \return each column is statedot vector
 */
Matd jointSpaceVelToStateDot(
    const Mat3& R_world, //Rotation matrix from body to world coords
    const Real roll,
    const Real pitch,
    const Matd& V)
{
  int rows = V.rows();
  int cols = V.cols();

  Matd S; //each column is a statedot vector
  S.resize(rows,cols); //same size as joint space velocity
  S.topRows(3) = R_world*V.block(3,0,3,cols); //rotate linear velocity into world coords
  S.block(3,0,3,cols) = velToEulerrateTransform(roll,pitch)*V.topRows(3); //transform angular vel (in body coords) to Euler angle rates
  S.bottomRows(rows-6) = V.bottomRows(rows-6); //copy joint rates

  return S;
}

void initWmrModelSurfaceContact(
    const WmrModel& mdl,
    const Surfaces& surfs,
    WmrState& state,
    std::vector<ContactFrame>& contacts,
    const OptimizationOptions& options)
{
  Vecd state_vec0 = state.toVecdEuler();

  Logical is_free(state_vec0.rows(), 0); //which elements of state vector are free to perturb
  is_free[2] = 1; //z
  is_free[3] = 1; //roll
  is_free[4] = 1; //pitch
  for (int frame_idx = 1; frame_idx < mdl.numFrames(); frame_idx++)
  {
    is_free[frame_idx + 5] = mdl.getFrame(frame_idx).is_fixed ? 0 : 1;
  }

  //get free elements of state vector
  Vecd x0;
  logicalIndexSrc(state_vec0, is_free, x0);

  Matd derr_dx; //derivative of error with respect to x.
  //computed in lambda function but required outside of lambda scope

  auto calcCostGradient = [&] (const Vecd& x, Real& cost, Vecd& gradient) mutable -> void
  {
    Vecd state_vec = state_vec0;
    std::vector<HTransform> HT_world;
    Vecd err;

    //calc cost
    {
      logicalIndexDst(x, is_free, state_vec);
      state.fromVecdEuler(state_vec);

      //calc HT_world
      std::vector<HTransform> HT_parent;
      stateToHTparent(mdl, state, HT_parent);
      HTparentToHTworld(mdl, HT_parent, HT_world);

      //calc contact geom
      calcWmrModelContactFrames(mdl, surfs, HT_world, contacts);

      //TODO, weight dz<0 error more than dz>0 error so that a subset of wheels can be in contact
      err.resize(contacts.size());
      for (size_t i = 0; i < contacts.size(); i++)
      {
        err[i] = contacts[i].dz;
      }
      //TODO, enforce additional holonomic joint constraints

      cost = err.dot(err);

      ROS_DEBUG_STREAM("x = " << x.transpose());
      ROS_DEBUG_STREAM("state = " << state_vec.transpose());
      ROS_DEBUG_STREAM("err = " << err.transpose());
    }

    //calc gradient
    {
      //calc Jacobian of contact point velocities with respect to joint space velocity
      Matd wheel_jacobians; //wheel jacobians
      calcWheelJacobians(mdl, HT_world, contacts, wheel_jacobians);

      //extract rows for z constraints
      Matd derr_djp; //derivative of error with respect to joint space perturbations
      derr_djp.resize(err.rows(),wheel_jacobians.cols());
      for (size_t i = 0; i < contacts.size(); i++)
      {
        derr_djp.row(i) = wheel_jacobians.row(i*3+2);
      }

      //convert to derivative of error with respect to state vector
      Real roll = state_vec[3];
      Real pitch = state_vec[4];
      Matd derr_dstate = jointSpaceVelToStateDot(
          HT_world[0].R, roll, pitch, derr_djp.transpose()).transpose();

      //copy the columns for free elements of state vector
      derr_dx.resize(err.rows(), x.rows());
      logicalIndexSrcCols(derr_dstate, is_free, derr_dx);

      ROS_DEBUG_STREAM("wheel_jacobians = \n" << wheel_jacobians);
      ROS_DEBUG_STREAM("derr_djp = \n" << derr_djp);
      ROS_DEBUG_STREAM("derr_dstate = \n" << derr_dstate);
      ROS_DEBUG_STREAM("derr_dx = \n" << derr_dx);

      //compute gradient
      gradient.resize(x.rows());
      gradient = 2.0*(err.transpose()*derr_dx).transpose();
    }
  };

  Vecd x = x0;
  Real cost;
  Vecd gradient;
  calcCostGradient(x, cost, gradient); //TODO, check if cost < cost_tol before computing gradient

  Real cost_prev = std::numeric_limits<Real>::infinity();
  for (int iter = 0; iter < options.max_iter; iter++)
  {
    ROS_DEBUG("iter = %d, cost = %f", iter, cost);
    ROS_DEBUG_STREAM("gradient = " << gradient.transpose());

    if (cost < options.cost_tol) { break; }
    if (cost_prev - cost < options.dcost_tol) { break; }

    //compute the Hessian
    Matd Hess = 2.0*derr_dx.transpose()*derr_dx;
    Vecd p = -Hess.llt().solve(gradient);

    ROS_DEBUG_STREAM("Hess = \n" << Hess);
    ROS_DEBUG_STREAM("p = " << p.transpose());

    Real dcost = gradient.dot(p); //expected decrease in cost
    if (std::abs(dcost) < options.dcost_tol) { break; } //gradient is near zero at minimum

    //TODO, linesearch to solve for alpha
    Real alpha = 0.9;
    x += alpha*p;
    calcCostGradient(x,cost,gradient);
  }
}

} //namespace
