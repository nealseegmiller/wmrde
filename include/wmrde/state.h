//state.h
//code for the state vector
//state vector comprises pose of the body frame (in world coords) and joint displacements
//state = [pose; joint displacements]
//pose = [orientation; position]
//orientation may be Euler angles or quaternion

//We also define joint space velocity and acceleration, qvel and qacc
//qvel is the same as the 1st derivative of state wrt time, except d/dt pose is replaced with the spatial velocity of the body frame (in body coords)
//qvel = [wx wy wz vx vy vz; joint rates]
//qacc is the same at the 2nd derivative of state wrt time, except d^2/dt^2 pose is replaced with spatial acceleration


#ifndef _WMRSIM_STATE_H_
#define _WMRSIM_STATE_H_

#include<WmrModel.h>

//for indexing state vector
#define SI_ORIENT 0
#define SI_POS SIZEORIENT
#define SI_JD (SIZEORIENT+3)

//VI_ANG and VI_LIN can't be swapped! spatial.h assumes this ordering
//these are for code readability only
#define VI_ANG 0
#define VI_LIN 3
#define VI_JR 6

#define NUMSTATE(nf) (nf - 1 + SIZEORIENT + 3)
#define NUMQVEL(nf) (nf - 1 + 6)

//convert frame index to state or qvel index, for fi > 0
#define TOSTATEI(fi) (fi - 1 + SI_JD)
#define TOQVELI(fi) (fi - 1 + VI_JR)

void stateToHT(const WmrModel &mdl, const Real state[], HomogeneousTransform HT_parent[], HomogeneousTransform HT_world[]);
void qvelToQdot(const int nf, const Real qvel[], const VecOrient orient, const Mat3 R_body_to_world, Real qdot[]);
void qdotToQvel(const int nf, const Real qdot[], const VecOrient orient, const Mat3 R_body_to_world, Real qvel[]);
void qvelToSpatialVel(const WmrModel &mdl, Mat6b Xup[], const Real qvel[], Vec6b v[]);

#endif //_WMRSIM_STATE_H_