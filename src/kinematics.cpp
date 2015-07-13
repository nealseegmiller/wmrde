#include <wmrde/kinematics.h>

//HT_world:		size nf, HT_world[i] is transform from frame i to world coords
//contacts:		size nw
//A:			3*npic x nv, Jacobian of contact point velocities with respect to joint space bias force
//return nc:	3*npic, number of rows
int wheelJacobians(const WmrModel& mdl, const HomogeneousTransform HT_world[], const WheelContactGeom contacts[], Real A[]) {


	//get from WmrModel
	const int nw = mdl.get_nw();; //number of wheels
	const int nf = mdl.get_nf(); //number of frames
	const int nv = NUMQVEL(nf); //size of joint space vel (cols of A)
	
	const int* wheelframeinds = mdl.get_wheelframeinds();
	const Frame* frames = mdl.get_frames();

	//
	int npic=0; //number of points in contact
	for (int wno = 0; wno < nw; wno++)
		npic = npic + contacts[wno].incontact;

	int nc = 3*npic; //number of constraints (rows of A)
	int row = 0;
	
	setMat(nc,nv,0.0,A); //init to zeros

	for (int wno = 0; wno < nw; wno++) { //wheel number
		if (contacts[wno].incontact) {

			int fi = wheelframeinds[wno]; //frame index
			Vec3 r; //translation vector

			while (fi > 0) {
				int dof_type = frames[fi].dof_type;
				int vi = TOQVELI(fi); //index in qvel

				if (dof_type < 3) { // 0,1,2 revolute
					//cross rotation axis with translation vector (from frame origin to contact pt)
					const Real* ax = HT_world[fi] + dof_type*SIZEVEC3; //rotation axis
					
					addmVec3(contacts[wno].HT_world+COL3, -1, HT_world[fi]+COL3, r);
					crossVec3(ax, r, A+S2I(row,vi,nc));
				} else { // 3,4,5 prismatic
					const Real* ax = HT_world[fi] + (dof_type-3)*SIZEVEC3;
					copyVec3(ax, A+S2I(row,vi,nc));
				}
				fi = frames[fi].parent_ind;
			}

			//body joint, fi=0
			addmVec3(contacts[wno].HT_world+COL3, -1, HT_world[0]+COL3, r);
			for (int ci=0; ci<3; ci++) { //column index
				const Real* ax = HT_world[0] + ci*SIZEVEC3; //axis
				crossVec3(ax, r, A + S2I(row,ci+VI_ANG,nc)); //angular
				copyVec3(ax, A + S2I(row,ci+VI_LIN,nc)); //linear
			}

			//rotate into contact frame coords
			for (int vi=0; vi<nv; vi++) {
				Real* A_ = A + S2I(row,vi,nc);
				multMatTVec3(contacts[wno].HT_world, A_, r);
				copyVec3(r, A_);
			}
			row += 3;
		}
	}

	return nc;
}


int trackJacobians(const WmrModel& mdl, const HomogeneousTransform HT_world[], const TrackContactGeom contacts[], Real A[]) {

	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//get from WmrModel
	const int nt = mdl.get_nt();
	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	const int* sprocketframeinds = mdl.get_sprocketframeinds();
	const Frame* frames = mdl.get_frames();
	
	//
	int npic = 0; //total number of points in contact
	for (int tno = 0; tno < nt; tno++) {
		npic += logicalCount(contacts[tno].get_np(), contacts[tno].incontact);
	}

	int nc = 3*npic; //number of constraints (rows of A)
	int row = 0; //row in A
	
	setMat(nc,nv,0.0,A); //init to zeros

	//for indexing spatial vectors
	const int ang = 0;
	const int lin = 3;

	for (int tno = 0; tno < nt; tno++) { //track number
		
		int incontactinds[ContactGeom::MAXNP];
		const int npic_track = logicalFind(contacts[tno].get_np(), contacts[tno].incontact, incontactinds);

		if (npic_track < 1)
			continue;

		int sprocket_fi = sprocketframeinds[tno];
		HomogeneousTransform HT_track_to_world;
		composeHT(HT_world[frames[sprocket_fi].parent_ind], frames[sprocket_fi].HT_parent_jd0, HT_track_to_world);

		Real rad = frames[sprocket_fi].rad;

		for (int i=0; i < npic_track; i++) {
			int pno = incontactinds[i]; //point number

			HomogeneousTransform HT_contact_to_world;
			composeHT(HT_track_to_world, contacts[tno].HT_track[pno], HT_contact_to_world);

			Vec3 r; //translation vector

			//sprocket angular rate
			mulcVec3(HT_contact_to_world + COL0, -rad, r); //temporary
			copyVec3(r, A + S2I(row,TOQVELI(sprocket_fi),nc));

			//traverse up the chain
			int fi = frames[sprocket_fi].parent_ind;

			while (fi > 0) {
				int dof_type = frames[fi].dof_type;
				int vi = TOQVELI(fi); //index in qvel

				if (dof_type < 3) { // 0,1,2 revolute
					//cross rotation axis with translation vector (from frame origin to contact pt)
					const Real* ax = HT_world[fi] + dof_type*SIZEVEC3; //rotation axis
					addmVec3(HT_contact_to_world+COL3, -1, HT_world[fi]+COL3, r);
					crossVec3(ax, r, A+S2I(row,vi,nc));
				} else { // 3,4,5 prismatic
					const Real* ax = HT_world[fi] + (dof_type-3)*SIZEVEC3; //translation axis
					copyVec3(ax, A+S2I(row,vi,nc));
				}
				fi = frames[fi].parent_ind;
			}

			//body joint, fi=0
			addmVec3(HT_contact_to_world+COL3, -1, HT_world[0]+COL3, r);
			for (int ci=0; ci<3; ci++) { //column index
				const Real* ax = HT_world[0] + ci*SIZEVEC3; //axis
				crossVec3(ax, r, A + S2I(row,ci+VI_ANG,nc)); //angular
				copyVec3(ax, A + S2I(row,ci+VI_LIN,nc)); //linear
			}

			
			//rotate into contact frame coords
			for (int vi=0; vi<nv; vi++) {
				Real* A_ = A + S2I(row,vi,nc);
				multMatTVec3(HT_contact_to_world, A_, r); //temporary
				copyVec3(r, A_);
			}

			row += 3;
		}
		
	}

	return nc;
}



//INPUTS
//mdl:		WmrModel object
//state:	size ns, state vector
//u:		size na, actuated joint rates
//HT_world: size nf, HT_world[i] transforms from frame i to world coords
//contacts: size nw or nt, ContactGeom objects
//OUTPUTS
//qvel:		size nv, joint space velocity 
//vc:		size nw, linear contact point velocities in contact coords
//TODO, handle tracks
void forwardVelKin(const WmrModel& mdl, const Real state[], const Real u[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, 
				   Real qvel[], Vec3 vc[]) {

	const int MAXNJ = WmrModel::MAXNF-1; //max number of joints
	const int MAXNV = NUMQVEL(WmrModel::MAXNF); //max size of qvel
	const int MAXNP = WmrModel::MAXNT * ContactGeom::MAXNP; //max number of contact points
	const int MAXNC = 3*MAXNP + WmrModel::MAXNJC; //max number of constraints


	const int nf = mdl.get_nf(); //number of frames in WmrModel obj
	const int nj = nf-1; //number of joints
	const int nv = NUMQVEL(nf); //size of qvel
	const Frame* frames = mdl.get_frames();

	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();

	const int* wheelframeinds = mdl.get_wheelframeinds();
	
	
	//contacts

	int np = 0; //total number of points
	bool incontact[MAXNP]; //concatenated, necessary?

	const WheelContactGeom* wcontacts =0;
	const TrackContactGeom* tcontacts =0;

	//to simplify indexing
	int whichtrack[MAXNP];
	Real dz[MAXNP];

	if (nw > 0) {
		wcontacts = static_cast<const WheelContactGeom*>(contacts);
		np = nw;
		for (int wno=0; wno < nw; wno++) {
			incontact[wno] = wcontacts[wno].incontact;
			dz[wno] = wcontacts[wno].dz;
		}
	} else if (nt > 0) {
		tcontacts = static_cast<const TrackContactGeom*>(contacts);
		for (int tno=0; tno < nt; tno++) {
			int np_ = tcontacts[tno].get_np(); //number of points for this track only
			copyVec(np_, tcontacts[tno].incontact, incontact+np);
			setVec(np_, tno, whichtrack+np);
			copyVec(np_, tcontacts[tno].dz, dz+np);
			np += np_;
		}
	}
	assert(np <= MAXNP); //DEBUGGING

	int incontactinds[MAXNP];
	int npic; //number of points in contact
	npic = logicalFind(np, incontact, incontactinds);


	//count constraints
	int ncc = 3*npic; //number of contact constraints
	int njc = mdl.get_njc(); //number of (holonomic) joint constraints
	int nc = ncc + njc; //total number of constraints

	//for indexing constraints (rows of A). indices of first:
	int contact_i0 = 0; //contact constraint
	int joint_i0 = ncc; //joint constraint

	//A*qvel = v;


	Real A[MAXNC * MAXNV]; //nc x nv
	Real v[MAXNC]; //nc x 1
	
	//contact constraints
	Real Acontact[3*MAXNP * MAXNV]; //ncc x nv
	setVec(ncc, 0.0, v+contact_i0); //init to zero

	if (nw > 0) {
		wheelJacobians(mdl, HT_world, wcontacts, Acontact);

		for (int i=0; i<npic; i++) {
			int wno = incontactinds[i];
			v[contact_i0 + (i*3) + 2] = -(wcontacts[wno].dz - mdl.dz_target[wno]) / mdl.tc_z[wno]; //Baumgarte's method
		}

	} else if (nt > 0) {
		trackJacobians(mdl, HT_world, tcontacts, Acontact);

		for (int i=0; i<npic; i++) {
			int pno = incontactinds[i]; //point number in list of all
			int tno = whichtrack[pno];

			v[contact_i0 + (i*3) + 2] = -(dz[pno] - mdl.dz_target[tno]) / mdl.tc_z[tno]; //Baumgarte's method
		}
	}

	copyMatBlock(ncc,0,0,ncc,nv,Acontact, nc,contact_i0,0,A);

	//holonomic joint constraints
	if (njc > 0) {
		Real c[WmrModel::MAXNJC];
		Real Jc[WmrModel::MAXNJC*MAXNJ];

		mdl.holonomicJointConstraints(mdl, state+SI_JD, 0, 
			c, Jc, 0, 0, 0);

		setMatBlock(nc,joint_i0,0,njc,VI_JR,0.0,A); //zeros for body DOF
		copyMatBlock(njc,0,0,njc,nj,Jc, nc,joint_i0,VI_JR,A);

		for (int i=0; i<njc; i++) //joint constraint number
			v[joint_i0 + i] = -c[i]/mdl.tc_j[i]; //Baumgarte's method
	}


	//DEBUGGING
	//std::cout << "A=\n"; printMatReal(nc,nv,A,-1,-1); std::cout << std::endl;
	//std::cout << "v=\n"; printMatReal(nc,1,v,-1,-1); std::cout << std::endl;


	//actuated (or fixed) joint rates
	int na = mdl.get_na();
	const int* actframeinds = mdl.get_actframeinds();

	bool vis_act[MAXNV]; //is fixed rate in joint space vel
	Real u_nv[MAXNV]; // size nv, padded

	setVec(nv,false,vis_act);
	
	for (int ai=0; ai < na; ai++) { //actuator index
		int vi = TOQVELI(actframeinds[ai]);
		vis_act[vi] = true;
		u_nv[vi] = u[ai];
	}

	//DEBUGGING
	
	//if wheel is not actuated and not in contact, make actuated to avoid rank deficiency
	for (int wno=0; wno < nw; wno++) {
		int fi = wheelframeinds[wno];
		if (!frames[fi].isactuated && !incontact[wno]) {
			int vi = TOQVELI(fi);
			vis_act[vi] = true;
			u_nv[vi] = 0.0; //set wheel vel = zero
			na++;
		}
	}
	

	//remove cols of A corresponding to fixed rates to the right hand side
	//Afree*qvel_free = b
	Real Afree[MAXNC * MAXNV]; // nc x (nv-na)
	Real qvel_free[MAXNV]; //(nv-na)
	Real b[MAXNC];
	

	//initialize
	copyVec(nc,v,b);

	int nfree=0;
	
	for (int vi = 0; vi < nv; vi++) {
		if ( vis_act[vi] ) {
			addmMatCol(nc,vi,A, 0,-u_nv[vi],b); //remove to right hand side
		} else {
			copyMatCol(nc,vi,A, nfree,Afree); //copy column
			nfree++;
		}
	}

	solve(nc,nv-na,Afree,b,qvel_free);

	nfree = 0;

	for (int vi=0; vi<nv; vi++) {
		if ( vis_act[vi] ) {
			qvel[vi] = u_nv[vi];
		} else {
			qvel[vi] = qvel_free[nfree];
			nfree++;
		}
	}

	//DEBUGGING
	//std::cout << "qvel=\n"; printMatReal(nv,1,qvel,-1,-1); std::cout << std::endl;

	if (vc != 0) {
		//output contact point velocities wrt ground
		Real vout[MAXNC];
		multMatVec(nc,nv,A,qvel,1.0, vout);
		int row = contact_i0;
		for (int pno=0; pno < np; pno++) {
			if (incontact[pno]) {
				copyVec3(vout+row, vc[pno]);
				row=row+3;
			} else {
				setcVec3(REALNAN, vc[pno]);
			}
		}

		//std::cout << "vc=\n"; printnVec3(nw,vc,-1,-1); //DEBUGGING
	}

}

//subfunction used by both odeKin and odeDyn
//use mdl.min_npic for kinematic sim, 0 for dynamic sim
void updateModelContactGeom(const WmrModel& mdl, const SurfaceVector& surfaces, const HomogeneousTransform HT_world[], const int min_npic, //inputs
	ContactGeom* contacts) { //outputs

	//get from WmrModel
	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	const Frame* frames = mdl.get_frames();

	//update contact geometry
	if (nw > 0) {
		WheelContactGeom* wcontacts = static_cast<WheelContactGeom*>(contacts);
		const int* wheelframeinds = mdl.get_wheelframeinds();

		for (int wno = 0; wno < nw; wno++) {
			int fi = wheelframeinds[wno];
			updateWheelContactGeom(surfaces, HT_world[fi], frames[fi].rad, wcontacts[wno]);
		}
		whichPointsInContactWheel(min_npic, nw, wcontacts);

	} else if (nt > 0) {
		TrackContactGeom* tcontacts = static_cast<TrackContactGeom*>(contacts);
		const int* sprocketframeinds = mdl.get_sprocketframeinds();

		for (int tno = 0; tno < nt; tno++) {
			int fi = sprocketframeinds[tno]; //frame index

			HomogeneousTransform HT_track_to_world;
			composeHT(HT_world[frames[fi].parent_ind], frames[fi].HT_parent_jd0, HT_track_to_world);

			updateTrackContactGeom(surfaces, HT_track_to_world, tcontacts[tno]);
		}
		whichPointsInContactTrack(min_npic, nt, tcontacts);
	}		
}

void odeKin(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, //inputs
	Real ydot[], HomogeneousTransform HT_parent[] ) { //outputs

	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//get from WmrModel
	const int nf = mdl.get_nf();

	//get control inputs
	Real u[WmrModel::MAXNA];
	mdl.controller(mdl, time, y, u, 0);

	//convert state to Homogeneous Transforms
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,y,HT_parent,HT_world);

	//update contact geometry
	updateModelContactGeom(mdl, surfaces, HT_world, mdl.min_npic, contacts);
	
	//compute joint space velocity
	Real qvel[MAXNV];
	forwardVelKin(mdl,y,u,HT_world,contacts,qvel,0);

	//convert to time derivative of state
	qvelToQdot(nf,qvel,y+SI_ORIENT,HT_world[0],ydot);

}



void initTerrainContact( const WmrModel mdl, const SurfaceVector& surfaces, ContactGeom* contacts, Real state[] ) {
	//TODO, handle quaternions
	assert(!WMRSIM_USE_QUATERNION);

	
	const int MAXNP = WmrModel::MAXNT * ContactGeom::MAXNP; //max number of contact points
	//max size of:
	const int MAXNS = NUMSTATE(WmrModel::MAXNF); //state vector
	const int MAXNE = MAXNP + WmrModel::MAXNJC; //error vector

	const int nf = mdl.get_nf();
	const int ns = NUMSTATE(nf);
	const int nw = mdl.get_nw(); //number of wheels
	const int nt = mdl.get_nt(); //number of tracks

	const Frame* frames = mdl.get_frames();

	//options
	const int max_iter = 15;
	const Real cost_tol = 1e-6;
	const Real dcost_tol = cost_tol/10;
	

	//which elements of state are free to change
	bool isfree[MAXNS]; // which elements of state vector are free
	
	setVec(ns,false,isfree);

	isfree[SI_ORIENT+0] = true; //rol
	isfree[SI_ORIENT+1] = true; //pit
	isfree[SI_POS+2] = true; //z
	for (int fi=1; fi < nf; fi++) 
		isfree[TOSTATEI(fi)] = !frames[fi].isfixed;

	//allocate vars needed outside of fCost, fGradient
	HomogeneousTransform HT_parent[WmrModel::MAXNF + WmrModel::MAXNW];
	HomogeneousTransform HT_world[WmrModel::MAXNF + WmrModel::MAXNW];
	
	//cast contacts
	WheelContactGeom* wcontacts;
	TrackContactGeom* tcontacts;
	if (nw > 0)
		wcontacts = static_cast<WheelContactGeom*>(contacts);
	else if (nt > 0)
		tcontacts = static_cast<TrackContactGeom*>(contacts);

	int ne; //size of error vector
	Real err[MAXNE];
	Real Jc[WmrModel::MAXNJC*(WmrModel::MAXNF-1)];
	Real derr_dx[MAXNE*MAXNS];

	Real p[MAXNS];
	Real x[MAXNS];
	Real cost;
	Real grad[MAXNS];

	//more
	int iter;
	
	Real cost_prev = REALMAX;
	Real Hess[MAXNS*MAXNS];
	Real HessL[MAXNS*MAXNS];

	//init x, comprises the free elements of state
	int nfree;
	nfree = logicalIndexIn( ns, isfree, state, x );

	//lambda closure, requires C++11
	auto initTerrainContactCost = [&] ( const Real x_[] ) mutable -> Real {
		//the following must be mutable: state, HT_parent, HT_world, ne, err, Jc
	
		logicalIndexOut( ns, isfree, x_, state );

		//convert state to Homogeneous Transforms
		stateToHT(mdl,state, HT_parent,HT_world);		

		updateModelContactGeom(mdl, surfaces, HT_world, mdl.min_npic, contacts);

		int npic = 0; //number of points in contact
		if (nw > 0) {
			for (int wno = 0; wno < nw; wno ++) {
				if (wcontacts[wno].incontact) {
					err[npic] = wcontacts[wno].dz - mdl.dz_target[wno];
					npic++;
				}
			}
		} else if (nt > 0) {
			for (int tno = 0; tno < nt; tno++) { //track number
				int np = tcontacts[tno].get_np(); //number of points for this track only
				for (int pno = 0; pno < np; pno++) {
					if (tcontacts[tno].incontact[pno]) {
						err[npic] = tcontacts[tno].dz[pno] - mdl.dz_target[tno];
						npic++;
					}
				}
			}
		}
		ne = npic;
		
		//additional (holonomic) joint constraints
		int njc = mdl.get_njc();
		if (njc > 0) {
			mdl.holonomicJointConstraints(mdl, state+SI_JD, 0, //inputs
				err+ne, Jc, 0, 0, 0); //outputs
			ne += njc;
		}

		assert(ne <= MAXNE);

		Real cost = dotVec(ne,err,err);

		return cost;

	};
	cost = initTerrainContactCost(x);

	//lambda closure, requires C++11
	auto initTerrainContactGradient = [&] ( Real grad_[] ) mutable {
		//must evaluate cost function before this
		//the following must be mutable: derr_dx

		//max size of:
		const int MAXNS = NUMSTATE(WmrModel::MAXNF); //state vector
		const int MAXNE = WmrModel::MAXNW + WmrModel::MAXNJC; //error vector
		const int MAXNV = NUMQVEL(WmrModel::MAXNF); //qvel vector

		Real A[3*WmrModel::MAXNW*MAXNV]; //contact Jacobian
		Real derrdot_dqvel[MAXNE*MAXNV]; //Jacobian of d/dt err wrt qvel
		Real derr_dstate[MAXNE*MAXNS]; //Jacobian of err wrt state

		int njc = mdl.get_njc();
		int npic = ne - njc; //number of points in contact
		int ncc = 3*npic; //number of contact constraints (rows of A)
		
		int nj = nf-1; //number of joints
		int nv = NUMQVEL(nf);

		if (nw > 0)
			wheelJacobians(mdl, HT_world, wcontacts, A);
		else if (nt > 0)
			trackJacobians(mdl, HT_world, tcontacts, A);

		//copy rows of A for z constraints to derrdot_dqvel
		for (int i=0; i<npic; i++) { //point in contact number
			int ri = i*3 + 2;
			copyMatRow(ncc,nv,ri,A, ne,i,derrdot_dqvel);
		}

		if (njc > 0) {
			//copy Jc to derrdot_dqvel
			setMatBlock(ne,npic,0,njc,VI_JR,0.0,derrdot_dqvel); //zeros for body DOF
			copyMatBlock(njc,0,0,njc,nj,Jc, ne,npic,VI_JR,derrdot_dqvel);
		}
	
		//convert to derr_dstate
		//TODO, eliminate temporary variables?
		Real qvel[MAXNV];
		Real qdot[MAXNS];
		for (int ri=0; ri<ne; ri++) { //row index
			copyRowToVec(ne,nv,ri,derrdot_dqvel, qvel);
			qvelToQdot(nf,qvel,state+SI_ORIENT,HT_world[0], qdot);
			copyVecToRow(qdot, ne,ns,ri,derr_dstate);
		}

		//copy columns from derr_dstate for free states
		int nfree=0;
		for (int si=0; si < ns; si++) { //state index
			if (isfree[si]) {
				copyMatCol(ne,si,derr_dstate, nfree,derr_dx);
				nfree++;
			}
		}

		//compute gradient
		for (int ci=0; ci < nfree; ci++) { //col index
			grad_[ci] = 2*dotVec(ne,err,derr_dx+S2I(0,ci,ne));
		}

		//DEBUGGING
		//std::cout << "A=\n"; printMatReal(ncc,nv,A,-1,-1); std::cout << std::endl;
		//std::cout << "Jc=\n"; printMatReal(njc,nj,Jc,-1,-1); std::cout << std::endl;
		//std::cout << "derrdot_dqvel=\n"; printMatReal(ne,nv,derrdot_dqvel,-1,-1); std::cout << std::endl;
		//std::cout << "derr_dstate=\n"; printMatReal(ne,ns,derr_dstate,-1,-1); std::cout << std::endl;
		//std::cout << "derr_dx=\n"; printMatReal(ne,nfree,derr_dx,-1,-1); std::cout << std::endl;
		//std::cout << "grad_=\n"; printMatReal(1,nfree,grad_,-1,-1); std::cout << std::endl;
	};


	for (iter=0; iter < max_iter; iter++) {
		
		
		//std::cout << "iter= " << iter << ", cost= " << cost << std::endl; //DEBUGGING
		//DEBUGGING, print contact height errors
		//if (nw > 0) {
		//	Real dz[WmrModel::MAXNW];
		//	for (int wno = 0; wno < nw; wno++)
		//		dz[wno] = wcontacts[wno].dz;
		//	std::cout << "dz = "; printMatReal(1,nw,dz,-1,-1);
		//}

		if (cost < cost_tol) 
			break;

		if (fabs(cost - cost_prev) < dcost_tol) 
			break;
		cost_prev = cost;

		if (iter==0)
			initTerrainContactGradient(grad);

		//compute Hessian, H = 2*(derr_dx^T * derr_dx)
		multMatTMat(ne, nfree, derr_dx, nfree, derr_dx, 2.0, Hess);

		chol(nfree,Hess,HessL);
		cholSolve(nfree,HessL,grad,p);

		mulcVec(nfree,-1.0,p); //negate

		Real gradp = dotVec(nfree,grad,p);
		if (fabs(gradp) < dcost_tol) {
			//gradient is zero at minimum, stop if close enough
			break;
		}

		Real alpha;
		alpha = linesearch( nfree, p, 1.0, initTerrainContactCost, initTerrainContactGradient, x, cost, grad);

		//DEBUGGING
		//std::cout << "grad=\n"; printMatReal(1,nfree,grad,-1,-1); std::cout << std::endl;
		//std::cout << "Hess=\n"; printMatReal(nfree,nfree,Hess,-1,-1); std::cout << std::endl;
		//std::cout << "p=\n"; printMatReal(nfree,1,p,-1,-1); std::cout << std::endl;
		//std::cout << "x=\n"; printMatReal(nfree,1,x,-1,-1); std::cout << std::endl;

	};


}
