#include <demo/talonmodel.h>

void talon(WmrModel& mdl, Real state[], Real qvel[]) {
	//dimensions from:
	//https://www.qinetiq-na.com/wp-content/uploads/brochure_talon.pdf

	const Real intom = 2.54/100;
	const Real lbtokg = .453592;

	//dimensions
	const Real rad = 5*intom; //sprocket radius
	const Real Wt = 6*intom; //width of track

	const int nt = 2;
	const Real L = (34*intom)-2*rad; //distance between forward/rear sprockets
	const Real B = (22.5*intom) - Wt; //distance between left/right track centers

	//masses
	const Real TotalMass = 115*lbtokg;
	const Real Mb = TotalMass; //mass of body
	const Real Ms = 1e-6; //mass of sprocket
	const Real Iy = .05*TotalMass*(rad*rad)/2; //rotational inertia of track about y axis

	//length, width, height of body
	const Real Lb = L;
	const Real Wb = B - Wt;
	const Real Hb = 1.5*rad;

	Vec3 cmb = {0,0,0}; //body center of mass (in body coords)
	Vec3 cms = {0,0,0}; //sprocket center of mass

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
	
	//left track
	setVec3(-L/2,B/2,0,pos);
	poseToHT(orient,pos,HT);
	mdl.addSprocketFrame("LeftTrack","Body",true,HT,rad,rad,L);

	//right track
	setVec3(-L/2,-B/2,0,pos);
	poseToHT(orient,pos,HT);
	mdl.addSprocketFrame("RightTrack","Body",true,HT,rad,rad,L);


	//set mass properties
	int nf = mdl.get_nf();
	const int* sprocketframeinds = mdl.get_sprocketframeinds();

	Mat3 I; //Inertia
	inertiaBox(Mb,Lb,Wb,Hb,I);
	mdl.setFrameMass(0,Mb,cmb,I);

	setcMat3(0.0,I);
	I[COL1+1] = Iy;
	for (int i=0; i<nt; i++)
		mdl.setFrameMass(sprocketframeinds[i],Ms,cms,I);


	//FOR KINEMATIC SIM
	mdl.min_npic = 3;

	//contact height error
	setVec(nt,-.01,mdl.dz_target);
	setVec(nt,.1,mdl.tc_z);

	//FOR DYNAMIC SIM

	//set function pointers
	mdl.controller = talonController;

	mdl.wheelGroundContactModel = uniformWgc;

	int npflat = 3;
	Real Kp = TotalMass*mdl.grav/(nt*npflat*-mdl.dz_target[0]);
	setWgcParams(Kp,mdl.wgc_p);

	mdl.actuatorModel = talonAct;
	mdl.act_p[0] = 2e3; //Kp
	mdl.act_p[1] = 0.0; //Ki
	mdl.act_p[2] = REALMAX; //max


	//ODE contact model parameters
	Real erp,cfm;
	KpKdToErpCfm(Kp, Kp/20, .04, erp, cfm);
	Real fds = 1.0/(Kp*1e-1);

	setVec(nt,erp,mdl.erp_z);
	setVec(nt,cfm,mdl.cfm_z);
	setVec(nt,fds,mdl.fds_x);
	setVec(nt,fds,mdl.fds_y);


	//initialize the state vector
	setEuler(DEGTORAD(0),DEGTORAD(0),DEGTORAD(0),euler);
	setVec3(0, -.9, rad - .01, pos);

#if WMRSIM_USE_QUATERNION
	eulerToQuat(euler,orient);
#else
	copyEuler(euler,orient);
#endif
	int ns = NUMSTATE(nf); //number of elements in state
	setVec(ns,0.0,state);
	copyOrient(orient,state+SI_ORIENT);
	copyVec3(pos,state+SI_POS);


	if (qvel != 0) { //not null
		//initialize qvel
		//setVec(mdl.get_na(), 0.0, qvel); //to zeros

		//to cmd
		Real u[WmrModel::MAXNA];
		mdl.controller(mdl, 0.0, state, u, qvel);
	}

}


void talonController(const WmrModel& mdl, const Real time, const Real state[], //inputs
					Real u[], Real qvel_cmd[]) { //outputs

	const Real intom = 2.54/100;

	//dimensions
	const Real rad = 5*intom; //sprocket radius
	const Real Wt = 6*intom; //width of track
	const Real B = (22.5*intom) - Wt; //distance between left/right track centers

	Real speed, omega; //commanded speed, yaw rate

	speed = 0.5;
	omega = 0;

	Real vl, vr; //vel left, vel right (m/s)

	vl = speed - (B/2)*omega;
	vr = speed + (B/2)*omega;

	u[0] = vl/rad;
	u[1] = vr/rad;

	
	if (qvel_cmd != 0) { //not null

		//get from WmrModel
		int nv = NUMQVEL(mdl.get_nf());
		const int nt = mdl.get_nt();
		const int* sprocketframeinds = mdl.get_sprocketframeinds();

		setVec(nv, 0.0, qvel_cmd);
		qvel_cmd[VI_ANG + 2] = omega;
		qvel_cmd[VI_LIN + 0] = speed;

		//sprocket velocities
		for (int i=0; i<nt; i++)
			qvel_cmd[TOQVELI(sprocketframeinds[i])] = u[i];
		
	}
}