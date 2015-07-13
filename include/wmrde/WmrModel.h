//WmrModel.h
//This class stores model information for Wheeled Mobile Robots



#ifndef _WMRSIM_WMRMODEL_H_
#define _WMRSIM_WMRMODEL_H_

//to write data to csv file
#include <iostream>
#include <fstream>

#include <assert.h>
#include <string>
#include <vector>


#include <algebra/spatial.h>

//Frame is a sub class of WmrModel
class Frame {
public:
	std::string name;
	int dof_type;		//degree of freedom type, 0-5
	int parent_ind;		//index of parent frame
	HomogeneousTransform HT_parent_jd0; //Homogeneous Transform of frame wrt parent if joint displacement=0

	bool iswheel;		//is a wheel frame
	bool issprocket;	//is a sprocket frame (for a track)
	//iswheel and issprocket cannot both be true
	bool isactuated;
	bool isfixed;		//if true, keep fixed when initializing terrain contact

	Real rad;			//radius of wheel or sprocket

	//for tracks
	Real rad2;
	Real L;

private:
	//mass properties
	Real mass;
	Vec3 cm;			//center of mass
	Mat3 I;				//moment of inertia

	Mat6b Is;			//spatial inertia
public:
	Frame();

	void setMass(const Real scalar_mass, const Vec3 center_of_mass, const Mat3 moment_of_inertia);
	static int convertDofTypeString(std::string dof_string); //converts string specifying joint type (revolute/prismatic) and axis (x,y,z) to integer

	//get methods
	const Real get_mass() const { return mass; }
	const Real* get_cm() const { return cm; }
	const Real* get_I() const { return I; }
	const Real* get_Is() const { return Is; }

};

inline int DofTypeToVec6bInd(const int dof_type) {
	//converts dof_type integer to index in Vec6b
	return (dof_type/3)*SIZEVEC3 + (dof_type % 3);
}

class WmrModel {

public:
	//buffer sizes
	static const int MAXNF = 13;		//max number of frames
	static const int MAXNW = 6;			//max number of wheels (or sprockets)
	static const int MAXNT = 2;			//max number of tracks (used to determine max number of contact points)
	static const int MAXNA = MAXNF-1;	//max number of actuated frames
	static const int MAXNJC = 1;		//max number of (holonomic) joint constraints
	static const int MAXNPAR = 20;		//max number of parameters
	

private:

	//SET AUTOMATICALLY BY ADDING FRAMES
	int nf;				//number of frames
	Frame frames[MAXNF];

	int njc;			//number of (holonomic) joint constraints
	
	
	//Dependent parameters, automatically updated
	int nw; //number of wheels
	int wheelframeinds[MAXNW]; //wheel frame index list
	int nt; //number of tracks
	int sprocketframeinds[MAXNW]; //sprocket frame index list
	int na; //number of actuated frames (or joints)
	int actframeinds[MAXNA]; //actuated frame index list
	
	Mat6b Is[MAXNF]; //spatial inertia

	void initFrame(const int fi);
	void update_wheelframeinds();
	void update_sprocketframeinds();
	void update_actframeinds();
	

public:
	//constructor
	WmrModel();

	//functions to construct model
	void addBodyFrame(const std::string name);
	void addJointFrame(const std::string name, const std::string parent_name, const std::string dof_string, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0);
	void addWheelFrame(const std::string name, const std::string parent_name, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0, 
		const Real radius);
	void addSprocketFrame(const std::string name, const std::string parent_name, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0, 
		const Real radius, const Real radius2, const Real len);
	
	void setFrameMass(const int fi, const Real mass, const Vec3 center_of_mass, const Mat3 moment_of_inertia);
	void set_isfixed(const int fi, const bool val) { frames[fi].isfixed = val; } 
	void set_njc( const int num_joint_constraints ) { njc = num_joint_constraints; }

	int nameToInd(const std::string Name) const; //returns index of frame by name

	//get methods
	//const correctness, const before { means read only access to WmrModel object
	const int get_nf() const { return nf; }
	
	
	const Frame* get_frames() const { return frames; }
	const Frame* get_frameptr(const int fi) const { return frames+fi; }
	//const Frame get_frame(const int fi) const { return frames[fi]; } //don't use this, slow!

	const int get_njc() const { return njc; }

	//
	const int get_nw() const { return nw; }
	const int* get_wheelframeinds() const { return wheelframeinds; }
	const int get_nt() const { return nt; }
	const int* get_sprocketframeinds() const { return sprocketframeinds; }
	const int get_na() const { return na; }
	const int* get_actframeinds() const { return actframeinds; }

	//TODO, make these private?
	//FOR KINEMATIC MODEL
	int min_npic;			//minimum number of points in contact
	Real dz_target[MAXNW];	//target Delta z (contact height error)
	//time constants for Baumgarte's stabilization method:
	Real tc_z[MAXNW];		//wheel-ground contact constraints, z dir
	Real tc_j[MAXNJC];		//holonomic joint constraints

	//FOR DYNAMIC MODEL
	Real grav;				//scalar acceleration of gravity
	bool use_constraints;


	//user must set these function pointers

	void (*controller) (const WmrModel& mdl, const Real time, const Real state[], //inputs
		Real u[], Real qvel_cmd[]); //outputs

	void (*wheelGroundContactModel) ( const int wheelno, const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
		Vec3 fw, Real J[]); //outputs
	Real wgc_p[MAXNPAR]; //wheel-ground contact model params
	//wheelno:	wheel number (in a size nw array)
	//J:		3x5 Jacobian of fw wrt [vx vy vz Rw dz]


	//one function for all actuators, may be coupled
	void (*actuatorModel) ( const Real params[], const Real ucmd[], const Real u[], const Real interr[], //inputs
		Real f[], Real err[], Real* dfdu); //outputs
	Real act_p[MAXNPAR]; //actuator model params
	//INPUTS
	//params:	size np
	//ucmd:		size na, commanded joint velocities
	//u:		size na, actual joint velocities
	//interr:	size na, integrated error
	//OUTPUTS
	//f:		size na, forces
	//err:		size na, error
	//dfdu:		na x na Jacobian, df/du

	//must also set njc > 0
	void (*holonomicJointConstraints) ( const WmrModel& mdl, const Real jd[], const Real jr[], //inputs
		Real c[], Real Jc[], Real f[], Real df_djd[], Real df_djr[]); //outputs
	//nj = number of joints
	//nc = number of joint constraints
	//INPUTS
	//jd:		size nj, joint displacement
	//jr:		size nj, joint rates (d/dt jd)
	//OUTPUTS
	//c:		size nc, evaluation of joint constraint function, zeros if satisfied
	//Jc:		nc x nj, Jacobian of c wrt jd
	//f:		size nc, force
	//df_djd:	nc x nj, Jacobian of f wrt jd
	//df_djr:	nc x nj, Jacobian of f wrt jr



	//the following are only used if no wheel-ground contact model is provided
	
	//for wheel-ground contact constraints
	Real erp_z[MAXNW];		//error reduction parameter, z dir (0-1)
	Real cfm_z[MAXNW];		//constraint force mixing, z dir (>0)
	Real fds_x[MAXNW];		//force dependent slip, x dir
	Real fds_y[MAXNW];		//force dependent slip, y dir

	//for holonomic joint constraints
	Real erp_j[MAXNJC];		//error reduction parameter
	Real cfm_j[MAXNJC];		//constraint force mixing

	
};

void inertiaBox(const Real m, const Real x, const Real y, const Real z, Mat3 I);
void inertiaCylinder(const Real m, const Real r, const Real h, const int axis, Mat3 I);
void KpKdToErpCfm( const Real kp, const Real kd, const Real dt, Real& erp, Real& cfm );

#endif //_WMRSIM_WMRMODEL_H_