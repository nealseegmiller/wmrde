#include <wmrde/wmrstate.h>
#include <wmrde/rosout.h>

namespace wmrde
{

inline Mat3 RotIdx(int idx, Real angle)
{
  if (idx == 0) { return Rotx(angle); }
  else if (idx == 1) { return Roty(angle); }
  else if (idx == 2 ){ return Rotz(angle); }
  else { return Mat3(); }
}

void stateToHTparent(
    const WmrModel& mdl,
    const WmrState& state,
    std::vector<HTransform>& HT_parent)
{
  assert(mdl.numFrames()-1 == state.joint_disp.rows());
  HT_parent.resize(mdl.numFrames());

  //get transform from body to world coords
  //TODO, check if use of .matrix() is correct
  HT_parent[0] = HTransform(state.orientation.matrix(), state.position);

  for (int frame_idx = 1; frame_idx < mdl.numFrames(); frame_idx++)
  {
    //abbreviate
    Real jd = state.joint_disp[frame_idx-1];
    int dof_type = mdl.getFrame(frame_idx).dof_type;
    const HTransform& HT_parent_jd0 = mdl.getFrame(frame_idx).HT_parent_jd0;

    if (dof_type < 3) //revolute joint
    {
      HT_parent[frame_idx].R = HT_parent_jd0.R * RotIdx(dof_type, jd);
      HT_parent[frame_idx].t = HT_parent_jd0.t; //copy translation
    }
    else //prismatic joint
    {
      HT_parent[frame_idx].R = HT_parent_jd0.R; //copy rotation
      Vec3 axis = HT_parent_jd0.R.row(dof_type-3); //translation axis
      HT_parent[frame_idx].t = HT_parent_jd0.t + jd*axis;
    }
  }
}

void HTparentToHTworld(
    const WmrModel& mdl,
    const std::vector<HTransform>& HT_parent,
    std::vector<HTransform>& HT_world)
{
  assert(mdl.numFrames() == HT_parent.size());

  HT_world.resize(HT_parent.size());
  HT_world[0] = HT_parent[0];

  for (int frame_idx = 1; frame_idx < mdl.numFrames(); frame_idx++)
  {
    int parent_idx = mdl.getFrame(frame_idx).parent_idx;
    HT_world[frame_idx] = HT_world[parent_idx].composeWith(HT_parent[frame_idx]);
  }
}

WmrState stepWmrState(
    const WmrState& state,
    const Vecd& joint_space_vel,
    const Real dt)
{
  assert(joint_space_vel.rows() == state.joint_disp.rows()+6);

  Mat3 R = state.orientation.matrix(); //rotation from body to world coords
  Vec6 body_vel = joint_space_vel.topRows(6); //spatial velocity of body frame
  int num_joints = state.joint_disp.rows();

  return WmrState(
      state.position + R*body_vel.bottomRows(3)*dt,
      stepRotationQuat(state.orientation, R*body_vel.topRows(3), dt),
      state.joint_disp + joint_space_vel.bottomRows(num_joints)*dt);
}

} //namespace

/*
//compute the spatial velocity of each frame from the joint space velocity
//using the RNEA, compare to jointSpaceBiasForce
//Xup:		size nf array of Plucker transforms, Xup[i] transforms spatial *motion* vector from parent(i) to i coords
//qvel:		size nv, joint space velocity
//v:		size nf, spatial velocity of each frame
//TODO, check this
void qvelToSpatialVel(const WmrModel &mdl, Mat6b Xup[], const Real qvel[], Vec6b v[]) {

	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	const Frame* frames = mdl.get_frames();

	//body frame
	copyArrayToVec6b(qvel,v[0]);

	for (int fi=1; fi < nf; fi++) {
		int parent_fi = frames[fi].parent_ind;
		int dof_type = frames[fi].dof_type;

		Vec6b vJ; //spatial joint velocity
		setcVec6b(0.0,vJ);
		vJ[DofTypeToVec6bInd(dof_type)] = qvel[TOQVELI(fi)];

		//v = Xup * v(parent) + vJ
		multPluckerVec6b(Xup[fi], v[parent_fi], v[fi]);
		addVec6b(v[fi],vJ,v[fi]);
	}
}
*/

