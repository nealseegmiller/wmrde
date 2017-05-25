#include <wmrde/kinematics.h>
#include <wmrde/util/index_util.h>
#include <wmrde/util/rosout.h>

namespace wmrde
{

void calcWheelJacobians(
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_world,
    const std::vector<ContactFrame>& contacts,
    Matd& wheel_jacobians)
{
  ROS_ASSERT(int(contacts.size()) == mdl.numWheelFrames()); //DEBUGGING

  //get number of points in contact
  int rows = contacts.size()*3;
  int cols = mdl.numDof();

  //allocate the wheel_jacobians matrix
  wheel_jacobians.resize(rows, cols);
  int row = 0;
  for (const ContactFrame& contact : contacts)
  {
    int frame_idx = contact.parent_idx; //index of parent wheel frame

    Eigen::Matrix<Real,3,Eigen::Dynamic> J;
    J.resize(3,cols);
    J.setZero();

    //first frame is wheel (or track) frame
    int vel_idx = frame_idx + 5; //index in joint space velocity
    Real radius = mdl.getFrame(frame_idx).wheel_geom->radius; //TODO, check for nullptr
    Vec3 direction = HT_world[frame_idx].R*contact.wheel_tangent_dir; //rotate into world coords
    J.col(vel_idx) = radius*direction; //rotate into world coords

    frame_idx = mdl.getFrame(frame_idx).parent_idx;

    while (frame_idx > 0)
    {
      vel_idx = frame_idx + 5; //index in joint space velocity
      int dof_type = mdl.getFrame(frame_idx).dof_type;
      Vec3 axis = HT_world[frame_idx].R.col(dof_type % 3);
      if (dof_type < 3) { //revolute joint
        //cross rotation axis with translation vector
        Vec3 r = contact.HT_world.t - HT_world[frame_idx].t; //from frame origin to contact pt
        J.col(vel_idx) = axis.cross(r);
      } else { //prismatic joint
        J.col(vel_idx) = axis;
      }
      frame_idx = mdl.getFrame(frame_idx).parent_idx;
    }

    //handle 6 Dof body joint (frame_idx == 0)
    Vec3 r = contact.HT_world.t - HT_world[frame_idx].t;
    for (vel_idx = 0; vel_idx < 3; vel_idx++)
    {
      Vec3 axis = HT_world[frame_idx].R.col(vel_idx);
      J.col(vel_idx) = axis.cross(r); //angular
      J.col(vel_idx+3) = axis; //linear
    }

    //rotate into contact frame coords and copy into output matrix
    wheel_jacobians.block(row,0,3,cols) = contact.HT_world.R.transpose() * J;
    row += 3;
  }
}

void forwardVelocityKinematics(
    const WmrModel& mdl,
    const WmrState& state,
    const Vecd& u,
    const std::vector<ContactFrame>& contacts,
    const std::vector<HTransform>& HT_world,
    JointSpaceVel& joint_space_vel)
{
  //count constraints
  int num_cc = 3*contacts.size();
  int num_hjc = 0; //TODO, mdl.numHolonomicJointConstraints();

  //Set up the constraint equation to solve for joint space velocity (qvel)
  //A*qvel = b

  Matd A(num_cc + num_hjc, mdl.numDof());
  Vecd b(A.rows());
  A.setZero();
  b.setZero();

  //contact constraints
  if (!contacts.empty())
  {
    Matd wheel_jacobians;
    calcWheelJacobians(mdl, HT_world, contacts, wheel_jacobians);

    //use Baumgarte's method to stabilize contact height error
    Vecd b_wheel(num_cc);
    b_wheel.setZero();
    for (size_t i = 0; i < contacts.size(); i++)
    {
      b_wheel[(i*3)+2] = (contacts[i].dz - mdl.dz_target)/mdl.dz_time_constant;
    }

    A.topRows(num_cc) = wheel_jacobians;
    b.topRows(num_cc) = b_wheel;
  }

  //TODO, holonomic joint constraints
  if (num_hjc > 0)
  {
    //TODO
//    mdl.holonomicJointConstraints(state.joint_disp,J,c); //TODO
//    A.bottomRows(num_hjc) = J;
//    b.bottomRows(num_hjc) = -c/mdl.hjc_time_constant;
  }
  ROS_DEBUG_STREAM("A = \n" << A);
  ROS_DEBUG_STREAM("b = \n" << b);

  //make logical for which elements of joint space velocity are free
  Logical is_free(mdl.numDof(), 1);
  for (int frame_idx : mdl.actuatedFrameIndices()) { is_free[frame_idx+5] = 0; }

  //TODO, if wheel is not actuated and not in contact, lock it to avoid rank deficiency

  Logical is_fixed = inverse(is_free);
  int num_free = countTrue(is_free);

  //make a padded u vector, the size of joint space velocity
  Vecd upad(mdl.numDof());
  upad.setZero();
  logicalIndexDst(u,is_fixed,upad);

  //remove cols of A corresponding to fixed rates to the right hand side
  //A_free * qvel_free = bp;
  Matd A_free(A.rows(), num_free);
  Vecd bp(b.rows()); //b + A_fixed * qvel_fixed

  int free_count = 0;
  for (int i = 0; i < mdl.numDof(); i++)
  {
    if (is_free[i]) {
      A_free.col(free_count++) = A.col(i); //copy column
    } else { //dof is fixed
      bp -= A.col(i)*upad[i]; //remove col to rhs
    }
  }

  Vecd qvel_free = A_free.householderQr().solve(bp);

  //copy to joint space velocity output
  joint_space_vel.resize(mdl.numDof());
  free_count = 0;
  for (int i = 0; i < mdl.numDof(); i++)
  {
    if (is_free[i]) { joint_space_vel[i] = qvel_free[free_count++]; }
    else { joint_space_vel[i] = upad[i]; }
  }
}

} //namespace

/*
void odeKin(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, //inputs
	Real ydot[], HomogeneousTransform HT_parent[] ) { //outputs

	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//get from WmrModel
	const int nf = mdl.get_nf();

	//get control inputs
	Real u[WmrModel::MAXNA];
	mdl.controller(mdl, time, y, u, 0);

	//convert state to Homogeneous Transforms
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,y,HT_parent,HT_world);

	//update contact geometry
	updateModelContactGeom(mdl, surfaces, HT_world, mdl.min_npic, contacts);
	
	//compute joint space velocity
	Real qvel[MAXNV];
	forwardVelKin(mdl,y,u,HT_world,contacts,qvel,0);

	//convert to time derivative of state
	qvelToQdot(nf,qvel,y+SI_ORIENT,HT_world[0],ydot);

}

*/
