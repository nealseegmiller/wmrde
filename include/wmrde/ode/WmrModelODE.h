//WmrModelODE.cpp
//This class stores an Open Dynamics Engine model of a WMR
//Used for benchmarking against WmrModel

#ifndef _WMRDE_WMRMODELODE_H_
#define _WMRDE_WMRMODELODE_H_

#include <assert.h>

#include <ode/ode.h>
#include <wmrde/algebra/transform.h>


//Open Dynamics Engine
class WmrModelODE {

public:
	//buffer sizes
	static const int MAXNB = 13; //max number of bodies
	static const int MAXNS = 10; //max number of surfaces
private:
	dWorldID world;
	dBodyID body[MAXNB];
	dJointID joint[MAXNB]; //joint[0] is null
	//use dJointGetBody(0) for parent
	//use dJointGetType(), returns dJointTypeHinge, dJointTypeBall
	dJointGroupID contactgroup;
	dContact contact;

	int nb; //number of bodies

	Vec3 pos_cm[MAXNB]; //pos_cm[i] position of center of mass in frame i coords (relative coordinates)
	//Open Dynamics Engine requires that center of mass coincide with point of reference


public:
	WmrModelODE(); //constructor
	~WmrModelODE(); //destructor

	void initParams();
	void addBodyBody(const Vec3 pos_center_of_mass); //add body for body frame
	void addBody(const int parent_ind, const int dof_type, const HomogeneousTransform HT_parent_jointdisp0, const Vec3 pos_center_of_mass);

	void setBodyMass( const int idx, const Real scalar_mass, const Mat3 moment_of_inertia);

	
	//get access
	const int get_nb() const { return nb; }
	const dWorldID get_world() const { return world; }
	const dBodyID* get_body() const { return body; }
	const dJointID* get_joint() const { return joint; }
	const dJointGroupID get_contactgroup() const { return contactgroup; }
	dSurfaceParameters* get_dSurfaceParameters() { return &contact.surface; }

	//GET POSE	
	//with respect to world frame
	void get_rot(const int idx, Mat3 R) const;
	void get_pos(const int idx, Vec3 pos) const;
	void get_HT(const int idx, HomogeneousTransform HT) const;
	
	Real get_jointdisp(const int idx) const;

	//SET POSE
	//to initialize simulation
	void set_rot(const int idx, const Mat3 R);
	void set_pos(const int idx, const Vec3 pos);
	void set_HT(const int idx, const HomogeneousTransform HT);


	//GET VEL
	//in world coords
	void get_angvel(const int idx, Vec3 angvel) const;
	void get_vel(const int idx, Vec3 vel) const;
	//in local coords
	void get_angvel_local(const int idx, Vec3 angvel) const;
	void get_vel_local(const int idx, Vec3 vel) const;
	
	Real get_jointrate(const int idx) const;

	//SET VEL
	//to initialize simulation
	//in world coords
	void set_angvel(const int idx, const Vec3 angvel);
	void set_vel(const int idx, const Vec3 vel);
	//in local coords
	void set_angvel_local( const int idx, const Vec3 angvel);
	void set_vel_local(const int idx, const Vec3 vel);
	
	void setMotorParams( const int idx, const Real MotorVel, const Real MotorFMax );
	void createContact( const int idx, const Real dz, const HomogeneousTransform HT_contact_to_world );


};


#endif

