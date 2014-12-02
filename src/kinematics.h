//kinematics.h
//functions for velocity kinematic motion prediction for WMRs
//and to initialize terrain contact
#ifndef _WMRSIM_KINEMATICS_H_
#define _WMRSIM_KINEMATICS_H_

#include <state.h>
//#include <contactgeom.h>
#include <collision.h> //necessary for initTerrainContact only
#include <eigensolve.h>


void wheelJacobians(const WmrModel& mdl, const HomogeneousTransform HT_world[], const WheelContactGeom contacts[], Real A[]);
void trackJacobians(const WmrModel& mdl, const HomogeneousTransform HT_world[], const TrackContactGeom contacts[], Real A[]);

void forwardVelKin(const WmrModel& mdl, const Real state[], const Real u[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, //inputs
				   Real qvel[], Vec3 vc[]); //outputs

void updateModelContactGeom(const WmrModel& mdl, const SurfaceVector& surfaces, const HomogeneousTransform HT_to_world[], const int min_npic, //inputs
	ContactGeom* contacts); //outputs

void odeKin(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, //inputs
	Real ydot[], HomogeneousTransform HT_to_parent[]); //outputs

void initTerrainContact( const WmrModel mdl, const SurfaceVector& surfaces, ContactGeom* contacts, Real state[] ); 

#endif //_WMRSIM_KINEMATICS_H_