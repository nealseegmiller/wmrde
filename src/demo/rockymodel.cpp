#include <demo/rockymodel.h>
//& to pass object by reference


//FUNCTIONS FOR THE ROCKY 7 ROVER
void rocky(WmrModel& mdl, Real state[], Real qvel[]) {
	//Rocky 7 rover model

	//dimensions
	const Real k1 = 10.5/100;		//vertical offset between (R)over reference to (D)ifferential
	const Real k2 = 12.075/100;	//forward offset between R and D
	const Real k3 = 20.0/100;		//horizontal distance between D and wheels
	const Real k4 = 28.8/100;		//distance from D to steering axis of front wheels
	const Real k5 = 12.5/100;		//height of D from wheel axles
	const Real k6 = 16.0/100;		//length of link from rocker joint to bogie joint
	const Real k7 = 6.9/100;		//length from bogie joint to front/rear bogie
	const Real k8 = 2.0/100;		//height of bogie joint from wheel axles
	const Real k9 = 139.0;		//angle of link between rocker and bogie joints (degrees)
	const Real k10 = 6.5/100;		//wheel radius

	const Real k6x = sin(DEGTORAD(k9-90.0))*k6;
	const Real k6z = cos(DEGTORAD(k9-90.0))*k6;
	

	const int nw=6; //number of wheels

	//masses (kg)
	const Real Mb = 11;		//mass of body
	const Real Mw = .5;		//mass of wheel
	const Real TotalMass = Mb + nw*Mw;

	const Real Lb = .40;		//length of body
	const Real Wb = .26;		//width of body
	const Real Hb = .16;		//height of body
	const Real rad = k10;
	const Real Ww = .08;		//wheel width

	Vec3 cmb = {18.0/100,0,15.5/100}; //body center of mass (in body coords)
	Vec3 cmw = {0,0,0}; //wheel center of mass

	//BUILD KINEMATIC TREE
	mdl.addBodyFrame("Body");

	HomogeneousTransform HT;
	VecEuler euler = {0,0,0};
	VecOrient orient;
#if WMRSIM_USE_QUATERNION
	eulerToQuat(euler,orient);
#else
	copyEuler(euler,orient);
#endif

	Vec3 pos;
	
	//LEFT SIDE
	//differential frame
	setVec3(k2,k3,k1,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("D1","Body","RY",false,HT);

	//steering frame
	setVec3(k4,0,0,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("S1","D1","RZ",true,HT);
	//wheel frame, front
	setVec3(0,0,-k5,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("A1","S1",true,HT,rad);

	//bogie frame
	setVec3(-k6x,0,-k6z,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("B1","D1","RY",false,HT);
	//wheel frame, bogie-front
	setVec3(k7,0,-k8,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("A3","B1",true,HT,rad);
	//wheel frame, bogie-rear
	setVec3(-k7,0,-k8,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("A5","B1",true,HT,rad);

	//RIGHT SIDE
	//differential frame
	setVec3(k2,-k3,k1,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("D2","Body","RY",false,HT);

	//steering frame
	setVec3(k4,0,0,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("S2","D2","RZ",true,HT);
	//wheel frame, front
	setVec3(0,0,-k5,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("A2","S2",true,HT,rad);

	//bogie frame
	setVec3(-k6x,0,-k6z,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("B2","D2","RY",false,HT);
	//wheel frame, bogie-front
	setVec3(k7,0,-k8,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("A4","B2",true,HT,rad);
	//wheel frame, bogie-rear
	setVec3(-k7,0,-k8,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("A6","B2",true,HT,rad);

	//set mass properties
	int nf = mdl.get_nf();
	const int* wheelframeinds = mdl.get_wheelframeinds();

	Mat3 I; //moment of inertia
	inertiaBox(Mb,Lb,Wb,Hb,I);
	mdl.setFrameMass(0,Mb,cmb,I);

	inertiaCylinder(Mw,rad,Ww,1,I);
	for (int i=0; i<nw; i++)
		mdl.setFrameMass(wheelframeinds[i],Mw,cmw,I);

	//FOR KINEMATIC MODEL

	mdl.min_npic = nw;

	setVec(nw,-.02,mdl.dz_target);
	setVec(nw,.1,mdl.tc_z);
	mdl.tc_j[0] = .1;

	//FOR DYNAMIC MODEL
	
	mdl.wheelGroundContactModel = uniformWgc;

	Real Kp = TotalMass*mdl.grav/(nw*-mdl.dz_target[0]);
	setWgcParams(Kp,mdl.wgc_p);

	mdl.actuatorModel = rockyAct;
	mdl.act_p[0] = 5; //Kp
	mdl.act_p[1] = 0; //Ki
	mdl.act_p[2] = REALMAX; //max

	//ODE contact model parameters
	Real erp,cfm;
	KpKdToErpCfm(Kp, Kp/20, .04, erp, cfm);
	Real fds = 1.0/(Kp*1e-1);

	setVec(nw,erp,mdl.erp_z);
	setVec(nw,cfm,mdl.cfm_z);
	setVec(nw,fds,mdl.fds_x);
	setVec(nw,fds,mdl.fds_y);
	mdl.erp_j[0] = .2;
	mdl.cfm_j[0] = 1e-6;

	//FOR BOTH
	//set function pointers
	mdl.controller = rockyController;
	mdl.holonomicJointConstraints = rockyConstraints;
	mdl.set_njc(1);

	//initialize the state vector
	setEuler(DEGTORAD(0),DEGTORAD(0),DEGTORAD(0),euler);
	setVec3(1.1, 0.0, (k10+k8+k6z-k1), pos);

#if WMRSIM_USE_QUATERNION
	eulerToQuat(euler,orient);
#else
	copyEuler(euler,orient);
#endif
	int ns = NUMSTATE(nf); //number of elements in state
	setVec(ns,0.0,state);
	copyOrient(orient,state+SI_ORIENT);
	copyVec3(pos,state+SI_POS);

	//initialize the joint space velocity
	if (qvel != 0) { //not null
		
		//setVec(mdl.get_na(), 0.0, qvel); //to zeros

		//to cmd
		Real u[WmrModel::MAXNA];
		mdl.controller(mdl, 0.0, state, u, qvel);
	}

}

void rockyController(const WmrModel& mdl, const Real time, const Real state[], //inputs
					Real u[], Real qvel_cmd[]) { //outputs

	//options

	//Dimensions
	//make sure these match
	Real k3 = 20.0/100;
	Real k4 = 28.8/100;
	Real k6 = 16.0/100;
	Real k9 = 139.0;
	Real k6x = sin(DEGTORAD(k9-90))*k6;
	Real k10 = 6.5/100;		//wheel radius

	Real L = k4+k6x; //length, front to rear
	Real w = k3; //half track width
	Real rad = k10; //wheel radius

	//indices of steer frames
	int S1_fi = 2;
	int S2_fi = 8;

	//for indexing u
	int S1_ai = 0;
	int A1_ai = 1;
	int A3_ai = 2;
	int A5_ai = 3;

	int S2_ai = 4;
	int A2_ai = 5;
	int A4_ai = 6;
	int A6_ai = 7;



	Real speed,turnrad; //commanded speed, turn radius

	speed = 0.5;
	turnrad = 1000;

	
	//if (time < 2.0) {
	//	speed = .5;
	//	turnrad = 1000;
	//} else if (time >= 2.0 && time < 6.0) {
	//	speed = .5;
	//	turnrad = -3.0;
	//} else {
	//	speed = .5;
	//	turnrad = 3.0;
	//}

	Real omega, gamma_l, gamma_r, vbl, vbr, vfl, vfr;

	if (fabs(turnrad) >= 1000) {
		omega = 0.0;

		gamma_l = 0.0;
		gamma_r = 0.0;

		vbl = speed;
		vbr = speed;
		vfl = speed;
		vfr = speed;

	} else {
		omega = speed/turnrad; //yaw rate

		//steer angles
		gamma_l = atan(L/(turnrad-w));
		gamma_r = atan(L/(turnrad+w));

		//steering limits
		Real lim = DEGTORAD(60);
		if (fabs(gamma_l) > lim)
			gamma_l = REALSIGN(gamma_l)*lim;
		if (fabs(gamma_r) > lim)
			gamma_r = REALSIGN(gamma_r)*lim;

		//wheel velocities
		//TODO, CHECK THIS!
		vbl = speed - omega*w; //back left
		vbr = speed + omega*w; //back right
		vfl = (speed - omega*w)*cos(gamma_l) + omega*L*sin(gamma_l); //front left
		vfr = (speed + omega*w)*cos(gamma_l) + omega*L*sin(gamma_l); //front right
	}

	//set u
	Real tc = .2;
	u[S1_ai] = 1/tc*(gamma_l - state[TOSTATEI(S1_fi)]);
	u[S2_ai] = 1/tc*(gamma_r - state[TOSTATEI(S2_fi)]);

	u[A1_ai] = vfl/rad;
	u[A3_ai] = vbl/rad;
	u[A5_ai] = vbl/rad;

	u[A2_ai] = vfr/rad;
	u[A4_ai] = vbr/rad;
	u[A6_ai] = vbr/rad;

	if (qvel_cmd != 0) { //not null

		//get from WmrModel
		const int nv = NUMQVEL(mdl.get_nf());
		const int na = mdl.get_na();
		const int* actframeinds = mdl.get_actframeinds();

		setVec(nv, 0.0, qvel_cmd);
		qvel_cmd[VI_ANG + 2] = omega;
		qvel_cmd[VI_LIN + 0] = speed;

		//wheel velocities
		for (int i=0; i<na; i++)
			qvel_cmd[TOQVELI(actframeinds[i])] = u[i];
	}
}


void rockyConstraints( const WmrModel& mdl, const Real jd[], const Real jr[], //inputs
	Real c[], Real Jc[], Real f[], Real df_djd[], Real df_djr[]) { //outputs

	//rocker joint indices
	const int D1_ji = 1-1;
	const int D2_ji = 7-1;

	c[0] = jd[D1_ji] + jd[D2_ji]; //angles should be equal magnitude, opposite sign
	
	
	const int nj = mdl.get_nf()-1; //number of joints
	

	//init to zeros
	setVec(nj,0.0,Jc);

	Jc[D1_ji] = 1.0;
	Jc[D2_ji] = 1.0;


	if ( f != 0 ) { //if not null
		//for dynamic sim
		Real Kp = 2250;
		Real Kd = Kp/20;
		Real cd = jr[D1_ji] + jr[D2_ji];

		f[0] = -(Kp*c[0] + Kd*cd);
		if ( df_djd != 0 ) {
			copyVec(nj,Jc,df_djd);
			mulcVec(nj,-Kp,df_djd);

			copyVec(nj,Jc,df_djr);
			mulcVec(nj,-Kd,df_djr);
		}
	}
}
