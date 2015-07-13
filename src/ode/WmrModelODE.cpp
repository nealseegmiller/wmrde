#include <wmrde/ode/WmrModelODE.h>

//constructor
WmrModelODE::WmrModelODE() {
	nb=0;

	dInitODE2(0);
	world = dWorldCreate();
	contactgroup = dJointGroupCreate (0);

}

WmrModelODE::~WmrModelODE() {

	dJointGroupDestroy (contactgroup);
	dWorldDestroy (world);
	dCloseODE();
}

void WmrModelODE::initParams() {

	//default parameter values

	//for world
	dWorldSetERP (world, 0.2); //default is 0.2, typical range is 0.1-0.8
	dWorldSetCFM (world, 1e-10); //default is 1e-10 for double, typical range is 10^-9 to 1
	dWorldSetGravity (world,0,0,-9.81);

	//for wheel-ground contact

	contact.surface.mode = dContactMu2 | dContactFDir1 | dContactSlip1 | dContactSlip2 | dContactApprox1;

	contact.surface.mu = dInfinity;
	contact.surface.mu2 = dInfinity;
	contact.surface.slip1 = 1e-4;
	contact.surface.slip2 = 1e-4;
	if (1) {
		contact.surface.mode = contact.surface.mode | dContactSoftERP | dContactSoftCFM;
		contact.surface.soft_erp = dWorldGetERP(world);
		contact.surface.soft_cfm = dWorldGetCFM(world);
	}
}

void WmrModelODE::addBodyBody( const Vec3 pos_center_of_mass ) {
	//add body for Wmr body

	assert(nb == 0);

	//add body frame
	int i = 0;
	nb++;

	copyVec3(pos_center_of_mass,pos_cm[i]); //must set this first!

	body[i] = dBodyCreate (world);
	
	HomogeneousTransform HT_world;
	setIdentityMat3(HT_world+0); //identity rotation
	setcVec3(0.0,HT_world+COL3); //at origin

	set_HT( i, HT_world );

	joint[i] = 0; //no joint for WMR body
}

void WmrModelODE::addBody( const int parent_ind, const int dof_type, const HomogeneousTransform HT_parent_jointdisp0, const Vec3 pos_center_of_mass) {
	//DofType 0-5 for [RX RY RZ PX PY PZ], Revolute or Prismatic about X,Y,Z axis

	assert(nb+1 <= MAXNB); //must not exceed allocated space
	assert(parent_ind < nb); //parent body must already be created

	int i = nb;
	nb++;

	copyVec3(pos_center_of_mass,pos_cm[i]); //must set this first!

	body[i] = dBodyCreate (world);

	//must set pose before creating joint, reference pose of parent body

	//HTParentJointDisp0 is HT_i_to_parent if joint displacement == 0;
	HomogeneousTransform HT_i_to_world;
	HomogeneousTransform HT_pi_to_world;

	get_HT(parent_ind,HT_pi_to_world);
	composeHT(HT_pi_to_world,HT_parent_jointdisp0,HT_i_to_world);

	set_HT(i,HT_i_to_world);

	//create and attach joint

	//joint axis
	int axisno = dof_type; //0-2 for X,Y,Z
	if (axisno > 3) axisno -= 3;
	Vec3 axis;
	setcVec3(0.0,axis);
	axis[axisno] = 1.0;

	
	if (dof_type >= 0 && dof_type <= 2) { //Revolute/Hinge joint

		joint[i] = dJointCreateHinge (world,0);
		
		dJointAttach (joint[i], body[i], body[parent_ind]); //order matters!

		Real* anchor = HT_i_to_world+COL3; //anchor at origin of frame i

		dJointSetHingeAnchor (joint[i], anchor[0], anchor[1], anchor[2]);
		dJointSetHingeAxis (joint[i], axis[0], axis[1], axis[2]);

	} else if (dof_type >= 3 && dof_type <= 5) { //Prismatic/Slider joint
		joint[i] = dJointCreateSlider (world,0);

		dJointAttach (joint[i], body[i], body[parent_ind]); //order matters!
		
		dJointSetSliderAxis (joint[i], axis[0], axis[1], axis[2]);
	}

}
void WmrModelODE::setBodyMass(const int idx, const Real scalar_mass, const Mat3 moment_of_inertia) {
	dMass M;
	const Real* I = moment_of_inertia;
	dMassSetParameters( &M, scalar_mass, 0, 0, 0, I[0], I[1+COL1], I[2+COL2], I[0+COL1], I[0+COL2], I[1+COL2]);
	dBodySetMass(body[idx], &M);

}

//GET POSE
void WmrModelODE::get_rot(const int idx, Mat3 R) const {
	//R:		rotation with respect to world frame

	const dReal* R_ = dBodyGetRotation(body[idx]);
	Mat3 RT;

	//copy, ODE's dMatrix3 has an extra row on bottom
	copyVec3(R_+0,RT+COL0);
	copyVec3(R_+4,RT+COL1);
	copyVec3(R_+8,RT+COL2);

	//convert row-major to column-major order
	copyTMat3(RT,R);
}

void WmrModelODE::get_pos(const int idx, Vec3 pos) const {
	//pos:		position of frame origin in world coords

	//account for center of mass offset
	dBodyGetRelPointPos(body[idx],-pos_cm[idx][0],-pos_cm[idx][1],-pos_cm[idx][2],pos);
}
void WmrModelODE::get_HT(const int idx, HomogeneousTransform HT) const {
	//HT:		Homogeneous Transform to world coordinates
	get_rot(idx, HT);
	get_pos(idx, HT+COL3);
}

Real WmrModelODE::get_jointdisp( const int idx ) const {
	int jt = dJointGetType(joint[idx]);
	Real val;
	if (jt == dJointTypeHinge) {
		val = dJointGetHingeAngle(joint[idx]);
	} else if (jt == dJointTypeSlider) {
		val = dJointGetSliderPosition(joint[idx]);
	}
	return val;
}

//SET POSE
void WmrModelODE::set_rot( const int idx, const Mat3 R ) {
	
	//convert column-major to row-major order
	Mat3 RT;
	copyTMat3(R,RT);
	
	dMatrix3 R_;
	
	//copy, ODE's dMatrix3 has an extra row on bottom
	copyVec3(RT+COL0,R_+0);
	copyVec3(RT+COL1,R_+4);
	copyVec3(RT+COL2,R_+8);
	
	dBodySetRotation(body[idx],R_);

}
void WmrModelODE::set_pos( const int idx, const Vec3 pos ) {
	
	dBodySetPosition (body[idx],pos[0],pos[1],pos[2]);

	//account center of mass offset
	//set orientation first!
	dVector3 pos_cm_world; //position of center of mass in world coords
	dBodyGetRelPointPos (body[idx],pos_cm[idx][0],pos_cm[idx][1],pos_cm[idx][2], pos_cm_world);
	dBodySetPosition (body[idx],pos_cm_world[0],pos_cm_world[1],pos_cm_world[2]);
}
void WmrModelODE::set_HT( const int idx, const HomogeneousTransform HT ) {
	
	//set orientation first to correctly account for center of mass offset
	set_rot(idx, HT); 
	set_pos(idx, HT+COL3);
}

//GET VEL
void WmrModelODE::get_angvel( const int idx, Vec3 angvel ) const {
	const dReal* ptr = dBodyGetAngularVel( body[idx] );
	copyVec3(ptr,angvel);
}
void WmrModelODE::get_vel( const int idx, Vec3 vel ) const {
	dBodyGetRelPointVel( body[idx], -pos_cm[idx][0], -pos_cm[idx][1], -pos_cm[idx][2], vel );
}
void WmrModelODE::get_angvel_local( const int idx, Vec3 angvel) const {
	Vec3 angvel_world;
	get_angvel(idx,angvel_world);

	//rotate to local coords
	Mat3 R_to_world;
	get_rot(idx,R_to_world);
	multMatTVec3(R_to_world,angvel_world,angvel);
}
void WmrModelODE::get_vel_local( const int idx, Vec3 vel) const {
	Vec3 vel_world;
	get_vel(idx,vel_world);

	//rotate to local coords
	Mat3 R_to_world;
	get_rot(idx,R_to_world);
	multMatTVec3(R_to_world,vel_world,vel);
}
Real WmrModelODE::get_jointrate( const int idx ) const {
	int jt = dJointGetType(joint[idx]);
	Real val;
	if (jt == dJointTypeHinge) {
		val = dJointGetHingeAngleRate(joint[idx]);
	} else if (jt == dJointTypeSlider) {
		val = dJointGetSliderPositionRate(joint[idx]);
	}
	return val;
}

//SET VEL
void WmrModelODE::set_angvel( const int idx, const Vec3 angvel ) {
	dBodySetAngularVel( body[idx], angvel[0], angvel[1], angvel[2] );
}
void WmrModelODE::set_vel( const int idx, const Vec3 vel ) {
	dBodySetLinearVel( body[idx], vel[0], vel[1], vel[2] );

	//account for center of mass offset
	//set angular velocity first!
	dVector3 vel_cm;
	dBodyGetRelPointVel( body[idx], pos_cm[idx][0], pos_cm[idx][1], pos_cm[idx][2], vel_cm );
	dBodySetLinearVel( body[idx], vel_cm[0], vel_cm[1], vel_cm[2] );
}
void WmrModelODE::set_angvel_local( const int idx, const Vec3 angvel ) {
	//rotate to world coords
	Vec3 angvel_world;
	Mat3 R_to_world;
	get_rot(idx,R_to_world);
	multMatVec3(R_to_world,angvel,angvel_world);

	set_angvel(idx,angvel_world);
}
void WmrModelODE::set_vel_local( const int idx, const Vec3 vel ) {
	//rotate to world coords
	Vec3 vel_world;
	Mat3 R_to_world;
	get_rot(idx,R_to_world);
	multMatVec3(R_to_world,vel,vel_world);

	set_vel(idx,vel_world);
}


void WmrModelODE::setMotorParams( const int idx, const Real MotorVel, const Real MotorFMax ) {
	dJointSetHingeParam (joint[idx], dParamVel, MotorVel);
	dJointSetHingeParam (joint[idx], dParamFMax, MotorFMax); //dInfinity, 0 turns off motor
}

void WmrModelODE::createContact( const int idx, const Real dz, const HomogeneousTransform HT_contact_to_world ) {
	//dz:		contact height error (< 0 indicates penetration)

	copyVec3(HT_contact_to_world+COL3,contact.geom.pos); //contact point = origin of contact frame
	copyVec3(HT_contact_to_world+COL2,contact.geom.normal); //normal = z axis
	
	contact.geom.depth = -dz; //depth should be > 0
	
	copyVec3(HT_contact_to_world+COL0,contact.fdir1); //friction direction 1 = x axis;
	
	dJointID cjoint = dJointCreateContact (world,contactgroup,&contact);
	//order matters!
	dJointAttach (cjoint,body[idx],0);
}


