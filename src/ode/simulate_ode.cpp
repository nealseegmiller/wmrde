#include <wmrde/ode/simulate_ode.h>

void convertToWmrModelODE( const WmrModel &mdl, WmrModelODE &mdl_ode ) {
	//convert WmrModel object to WmrModelODE object
	const int nf = mdl.get_nf();
	const Frame* frames = mdl.get_frames();

	mdl_ode.initParams();

	//kinematics
	mdl_ode.addBodyBody( frames[0].get_cm() ); //body frame
	for (int fi = 1; fi < nf; fi++) {
		mdl_ode.addBody( frames[fi].parent_ind, frames[fi].dof_type, frames[fi].HT_parent_jd0, frames[fi].get_cm() );
	}

	//mass properties
	//set mass
	for (int fi = 0; fi < nf; fi++) {
		if (frames[fi].get_mass() > 0)
			mdl_ode.setBodyMass( fi, frames[fi].get_mass(), frames[fi].get_I() );
	}

	//contact properties?
}

void getHTODE( const WmrModel &mdl, const WmrModelODE &mdl_ode, HomogeneousTransform HT_world[], HomogeneousTransform HT_parent[] ) {
	const int nf = mdl.get_nf();
	

	//get HT wrt world
	for (int fi = 0; fi < nf; fi++)
		mdl_ode.get_HT(fi,HT_world[fi]);

	if (HT_parent != 0) {
		//get HT wrt parent
		const Frame* frames = mdl.get_frames();

		copyHT(HT_world[0],HT_parent[0]); //body frame
		for (int fi = 1; fi < nf; fi++) {
			int parent_fi = frames[fi].parent_ind;
			composeInvHT(HT_world[parent_fi],HT_world[fi],HT_parent[fi]);
		}
	}
}

void getStateODE( const WmrModelODE &mdl_ode, Real state[] ) {
	//pose of the WMR body
	HomogeneousTransform HT_body_to_world;
	mdl_ode.get_HT(0,HT_body_to_world);
	HTToPose(HT_body_to_world,state+SI_ORIENT,state+SI_POS);

	//joint displacements
	const int nb = mdl_ode.get_nb();
	for (int i=1; i<nb; i++)
		state[SI_JD+i-1] = mdl_ode.get_jointdisp(i);
}


void setStateODE( const WmrModel &mdl, const Real state[], WmrModelODE& mdl_ode ) {
	//assumes kinematics of WmrModel and WmrModelODE objects are consistent!

	HomogeneousTransform HT_parent[WmrModelODE::MAXNB];
	HomogeneousTransform HT_world[WmrModelODE::MAXNB];

	stateToHT( mdl, state, HT_parent, HT_world );

	const int nb = mdl_ode.get_nb();
	for (int i=0; i<nb; i++)
		mdl_ode.set_HT(i,HT_world[i]);

}

//joint space velocity
void getQvelODE( const WmrModelODE &mdl_ode, Real qvel[] ) {
	//spatial velocity of body frame
	mdl_ode.get_angvel_local(0,qvel+VI_ANG);
	mdl_ode.get_vel_local(0,qvel+VI_LIN);

	//joint rates
	const int nb = mdl_ode.get_nb();
	for (int i=1; i<nb; i++)
		qvel[VI_JR+i-1] = mdl_ode.get_jointrate(i);
}

void setQvelODE( const WmrModel &mdl, const Real qvel[], WmrModelODE &mdl_ode ) {
	const int nf = mdl.get_nf();
	
	//get Plucker transforms, Xup
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	getHTODE( mdl, mdl_ode, HT_world, HT_parent );

	Mat6b Xup[WmrModel::MAXNF];
	for (int fi=0; fi<nf; fi++)
		invHTToPlucker(HT_parent[fi],Xup[fi]);

	Vec6b vframe[WmrModel::MAXNF]; //spatial velocity of each frame, local coords
	qvelToSpatialVel(mdl, Xup, qvel, vframe);

	for (int fi=0; fi<nf; fi++) {		
		mdl_ode.set_angvel_local(fi, vframe[fi]+0); //set angular velocity first!
		mdl_ode.set_vel_local(fi, vframe[fi]+COL1);
	}
}

void stepODE( const Real time, const WmrModel& mdl, const SurfaceVector& surfaces, const Real dt, //inputs
	WmrModelODE& mdl_ode, HomogeneousTransform HT_parent[], WheelContactGeom contacts[]) { //outputs

	//TODO
	//compare to odeDyn()

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nw = mdl.get_nw();
	const int na = mdl.get_na();

	const int* wheelframeinds = mdl.get_wheelframeinds();
	const int* actframeinds = mdl.get_actframeinds();

	const Frame* frames = mdl.get_frames();


	HomogeneousTransform HT_world[WmrModel::MAXNF];

	Real state[MAXNS];
	getStateODE(mdl_ode, state);

	//stateToHT(mdl, state, HT_parent, HT_world); //DEBUGGING

	//get HT wrt world, parent
	getHTODE( mdl, mdl_ode, HT_world, HT_parent );

	//wheel-ground contact joints

	
	dSurfaceParameters* sp = mdl_ode.get_dSurfaceParameters();

	for (int wno = 0; wno < nw; wno++) {
		int fi = wheelframeinds[wno];
		updateWheelContactGeom(surfaces, HT_world[fi], frames[fi].rad, contacts[wno] );

		//set surface parameters according to WmrModel?
		sp->soft_erp = mdl.erp_z[wno];
		sp->soft_cfm = mdl.cfm_z[wno];
		sp->slip1 = mdl.fds_x[wno];
		sp->slip2 = mdl.fds_y[wno];

		if (contacts[wno].dz <= 0)
			mdl_ode.createContact( fi, contacts[wno].dz, contacts[wno].HT_world );
	}

	//set motor parameters
	
	//get control inputs
	Real ucmd[WmrModel::MAXNA];
	mdl.controller(mdl, time, state, ucmd, 0);
	
	for (int ai=0; ai<na; ai++) {
		int fi = actframeinds[ai];
		mdl_ode.setMotorParams(fi,ucmd[ai],dInfinity);
	}
	
	dWorldStep (mdl_ode.get_world(),dt);
	dJointGroupEmpty (mdl_ode.get_contactgroup());

}

