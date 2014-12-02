//dynamics.h
//functions for dynamic motion prediction for WMRs

#ifndef _WMRSIM_DYNAMICS_H_
#define _WMRSIM_DYNAMICS_H_

#include <kinematics.h>

void subtreeInertias(const WmrModel &mdl, const Mat6b Xup[], Mat6b Is_subt[]);
void jointSpaceInertia( const WmrModel& mdl, const Mat6b Xup[], const Mat6b Is_subt[], Real H[] );
void jointSpaceBiasForce(const WmrModel& mdl, const Mat6b Xup[], const Real qvel[], Real C[]);

void forwardDyn(const WmrModel& mdl, const Real state0[], const Real qvel0[], const Real u_cmd[], const Real u_interr[], 
	const HomogeneousTransform HT_parent[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, const Real dt, //inputs
	Real qacc[], Real u_err[]); //outputs

void odeDyn(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, const Real dt, //inputs
	Real ydot[], HomogeneousTransform HT_to_parent[]); //outputs

#endif //_WMRSIM_DYNAMICS_H_