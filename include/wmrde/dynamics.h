//dynamics.h
//functions for dynamic motion prediction for WMRs

#ifndef _WMRDE_DYNAMICS_H_
#define _WMRDE_DYNAMICS_H_

#include <wmrde/kinematics.h>

//classes to reduce the number of function parameters
//always pass by reference!

//if MAXNP large, requires large stack size (10 Mb)
//http://msdn.microsoft.com/en-US/library/8cxs58a6(v=vs.120).aspx
//Configuration Properties -> Linker -> System -> Stack Reserve Size
//TODO, allocate large matrices on heap?

class ConstraintJacobian {
public:
	//for allocation
	static const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	static const int MAXNP = WmrModel::MAXNT * ContactGeom::MAXNP; //max number of contact points
	static const int MAXNC = 3*MAXNP + WmrModel::MAXNA + WmrModel::MAXNJC;

	//
	Real contact[3*MAXNP * MAXNV];
	Real act[WmrModel::MAXNA * MAXNV]; //actuator
	Real joint[WmrModel::MAXNJC * MAXNV]; //holonomic joint constraints
	Real all[MAXNC * MAXNV];
};

class ControllerIO {
public:
	//input
	Real cmd[WmrModel::MAXNA]; //commanded
	Real interr[WmrModel::MAXNA]; //integrated error
	//output
	Real err[WmrModel::MAXNA]; //error: commanded - actual
};


void subtreeInertias(const WmrModel &mdl, const Mat6b Xup[], Mat6b Is_subt[]);
void jointSpaceInertia( const WmrModel& mdl, const Mat6b Xup[], const Mat6b Is_subt[], Real H[] );
void jointSpaceBiasForce(const WmrModel& mdl, const Mat6b Xup[], const Real qvel[], Real C[]);

void HandC(const WmrModel& mdl, const HomogeneousTransform HT_parent[], const Real qvel[], //input
	Real H[], Real C[]); //output

void forwardDyn(const WmrModel& mdl, const Real state0[], const Real qvel0[], ControllerIO& u, 
	const HomogeneousTransform HT_parent[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, const Real dt, //inputs
	Real qacc[]); //outputs

//TODO change name?
void constraintJacobians(const WmrModel& mdl, const Real state0[], const Real qvel0[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, //input
	ConstraintJacobian& A ); //output

int parseWheelContacts(const WheelContactGeom* wcontacts, const int nw, bool incontact[], Real dz0[]);
int parseTrackContacts(const TrackContactGeom* tcontacts, const int nt,	bool incontact[], Real dz0[], int whichtrack[]);

//TODO, remove this
void forwardDynErpCfm(const WmrModel& mdl, const Real state0[], const Real qvel0[], const Real u_cmd[], const ContactGeom* contacts, const Real dt, 
	Real H[], const Real C[], const ConstraintJacobian& A, //input
	Real qacc[]); //output

void forwardDynForceBalance(const WmrModel& mdl, const Real state0[], const Real qvel0[], ControllerIO& u, const ContactGeom* contacts, const Real dt, 
	const Real H[], const Real C[], const ConstraintJacobian& A, //input
	Real qacc[]); //output

//TODO forwardDynUnc()

void odeDyn(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, const Real dt, //inputs
	Real ydot[], HomogeneousTransform HT_parent[]); //outputs

#endif //_WMRDE_DYNAMICS_H_
