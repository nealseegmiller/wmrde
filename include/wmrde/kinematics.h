//kinematics.h
//functions for kinematic motion prediction for WMRs
//and to initialize terrain contact

#ifndef _WMRDE_KINEMATICS_H_
#define _WMRDE_KINEMATICS_H_

#include <wmrde/state.h>
#include <wmrde/collision.h> //necessary for initTerrainContact only
#include <wmrde/eigensolve.h>


int wheelJacobians(const WmrModel& mdl, const HomogeneousTransform HT_world[], const WheelContactGeom contacts[], Real A[]);
int trackJacobians(const WmrModel& mdl, const HomogeneousTransform HT_world[], const TrackContactGeom contacts[], Real A[]);

void forwardVelKin(const WmrModel& mdl, const Real state[], const Real u[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, //inputs
				   Real qvel[], Vec3 vc[]); //outputs

void updateModelContactGeom(const WmrModel& mdl, const SurfaceVector& surfaces, const HomogeneousTransform HT_world[], const int min_npic, //inputs
	ContactGeom* contacts); //outputs

void odeKin(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, //inputs
	Real ydot[], HomogeneousTransform HT_parent[]); //outputs

void initTerrainContact( const WmrModel mdl, const SurfaceVector& surfaces, ContactGeom* contacts, Real state[] ); 

#endif //_WMRDE_KINEMATICS_H_
