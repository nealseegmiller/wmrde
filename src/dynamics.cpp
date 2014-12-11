#include <dynamics.h>

//Xup:		size nf array of Plucker transforms, Xup[i]^T*f transforms spatial *force* vector f from frame i to parent coords
//Is_subt:	size nf array of spatial inertias, Is_subt[i] is the inertia of the subtree rooted at frame i
void subtreeInertias(const WmrModel &mdl, const Mat6b Xup[], Mat6b Is_subt[]) {

	
	const int nf = mdl.get_nf();
	const Frame* frames = mdl.get_frames();

	//init to frame spatial inertia
	for (int fi=0; fi < nf; fi++) { 
		copyMat6b(frames[fi].get_Is(),Is_subt[fi]);
	}
	
	Mat6b Is_pc;

	for (int fi=(nf-1); fi > 0; fi--) { //frame index
		//transform to parent coords

		//multPluckerTMat6bPlucker(Xup[fi], Is_subt[fi], Is_pc); 
		multPluckerTInertiaPlucker(Xup[fi], Is_subt[fi], Is_pc); 
		
		int parent_fi = frames[fi].parent_ind; //parent frame index
		addMat6b(Is_subt[parent_fi], Is_pc, Is_subt[parent_fi]); //add to parent inertia
	}
	
}

//compute joint space inertia using Composite Rigid Body Algorithm
//adapted from HandC.m in spatial_v2 library
//http://users.cecs.anu.edu.au/~roy/spatial/
void jointSpaceInertia( const WmrModel& mdl, const Mat6b Xup[], const Mat6b Is_subt[], Real H[] ) {


	//from WmrModel
	const int nf = mdl.get_nf();
	const int nv = TOQVELI(nf);
	const Frame* frames = mdl.get_frames();


	setMat(nv,nv,0.0,H); //init to zeros

	//copy Is_subt for body frame to 6x6 block of H
	Real H_[6*6];
	copyMat6bToArray(Is_subt[0],H_);
	copyMatBlock(6,0,0,6,6,H_, nv,0,0,H);


	Vec6b force_h; //spatial force vector
	Vec6b tmp;
	int dof_type;

	for (int fi=1; fi<nf; fi++) { //frame index
		int vi = TOQVELI(fi); //vel index (in qvel)
		dof_type = frames[fi].dof_type;

		copyMat6bColToVec6b(dof_type,Is_subt[fi],force_h);
		H[S2I(vi,vi,nv)] = force_h[dof_type]; //diagonal term

		int fj = fi; //2nd frame index
		//move up the chain to compute off-diagonal terms
		while (fj > 0) {
			multPluckerTVec6b(Xup[fj],force_h,tmp);
			copyVec6b(tmp,force_h);
			fj = frames[fj].parent_ind;

			if (fj == 0) {
				//body joint
				copyVec6bToArray(force_h, H+S2I(0,vi,nv));
				copyVecToRow(H+S2I(0,vi,nv), nv,6,vi,H); //symmetry
			} else {
				int vj = TOQVELI(fj); //2nd vel index
				dof_type = frames[fj].dof_type;

				H[S2I(vj,vi,nv)] = force_h[DofTypeToVec6bInd(dof_type)];
				H[S2I(vi,vj,nv)] = H[S2I(vj,vi,nv)]; //symmetry
			}
		}
	}
}

//compute joint space bias force using Recursive Newton-Euler Algorithm
//adapted from HandC.m in spatial_v2 library
//http://users.cecs.anu.edu.au/~roy/spatial/
void jointSpaceBiasForce(const WmrModel& mdl, const Mat6b Xup[], const Real qvel[], Real C[]) {

	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	const Frame* frames = mdl.get_frames();

	//variable names from HandC.m in spatial_v2 library
	//spatial vectors for each frame:
	Vec6b v[WmrModel::MAXNF]; //velocity
	Vec6b avp[WmrModel::MAXNF]; //accleration (includes Coriolis, centripetal, and gravity terms)
	Vec6b fvp[WmrModel::MAXNF]; //force

	//body frame
	copyArrayToVec6b(qvel,v[0]);
	//spatial accel of gravity in body coords
	Vec6b a_grav;
	setcVec3(0.0, a_grav); //angular
	mulcVec3(Xup[0]+COL2, mdl.grav, a_grav+COL1); //linear, world z axis in body coords * grav acc

	//accelerating the base is faster than explicitly adding gravity force to each body
	copyVec6b(a_grav, avp[0]);
	
	Vec6b tmp;

	for (int fi=0; fi < nf; fi++) {
		if (fi > 0) {
			int parent_fi = frames[fi].parent_ind;
			int dof_type = frames[fi].dof_type;

			
			Vec6b vJ; //spatial joint velocity
			setcVec6b(0.0,vJ);
			vJ[DofTypeToVec6bInd(dof_type)] = qvel[TOQVELI(fi)];

			//v = Xup * v(parent) + vJ
			multPluckerVec6b(Xup[fi], v[parent_fi], v[fi]);
			addVec6b(v[fi],vJ,v[fi]);

			//avp = Xup * avp(parent) + v x vJ
			multPluckerVec6b(Xup[fi], avp[parent_fi], avp[fi]);
			crossVec6bMotion(v[fi],vJ,tmp);
			addVec6b(avp[fi],tmp,avp[fi]);
		}
		if (frames[fi].get_mass() > 0) {
			//fvp = Is*avp + v x Is*v
			const Real* Is = frames[fi].get_Is();

			multMatVec6b(Is,v[fi],tmp);
			crossVec6bForce(v[fi],tmp,fvp[fi]);
			multMatVec6b(Is,avp[fi],tmp);
			addVec6b(fvp[fi],tmp,fvp[fi]);
		} else {
			setcVec6b(0.0, fvp[fi]);
		}
	}

	for (int fi=nf-1; fi>0; fi--) {
		int parent_fi = frames[fi].parent_ind;
		int dof_type = frames[fi].dof_type;

		C[TOQVELI(fi)] = fvp[fi][DofTypeToVec6bInd(dof_type)];

		multPluckerTVec6b(Xup[fi],fvp[fi],tmp);
		addVec6b(fvp[parent_fi],tmp,fvp[parent_fi]);
	}
	copyVec6bToArray(fvp[0],C);

	//DEBUGGING
	//std::cout << "v=\n"; printnVec6b(nf,v,-1,-1);
	//std::cout << "avp=\n"; printnVec6b(nf,avp,-1,-1);
	//std::cout << "fvp=\n"; printnVec6b(nf,fvp,-1,-1);
	//std::cout << "C=\n"; printMatReal(nv,1,C,-1,-1);

	return;
}


//compute joint space inertia and bias force
void HandC(const WmrModel& mdl, const HomogeneousTransform HT_parent[], const Real qvel[], //input
	Real H[], Real HL[], Real C[]) { //output

	//OPTIONS
	const bool do_reuse_H = false; //reuse joint space inertia
	const bool do_approx_C = false; //treat WMR as single rigid body when computing joint space bias force

	//for allocation
	const int MAXNV = NUMQVEL(WmrModel::MAXNF); //max size of qvel

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	//static vars
	//TODO, eliminate static vars, store in WmrModel instead
	static Mat6b Is_subt[WmrModel::MAXNF];
	static Real H_prev[MAXNV*MAXNV]; 
	static Real HL_prev[MAXNV*MAXNV]; //lower triangular
	static bool H_init = false; //flag, H is initialized

	//convert HomogeneousTransforms to Plucker transforms
	//Xup[i] transforms spatial *motion* vector from frame parent(i) to i coords
	//Xup[i]^T transforms spatial *force* vector from frame i to parent(i) coords
	Mat6b Xup[WmrModel::MAXNF];
	for (int fi=0; fi<nf; fi++) 
		invHTToPlucker(HT_parent[fi],Xup[fi]);

	//joint space inertia, H
	
	//Is_subt[i] is the 6x6 spatial inertia of the subtree rooted at frame i, in frame i coords.
	
	if (!do_reuse_H || !H_init) {
		subtreeInertias(mdl, Xup, Is_subt); //expensive!
		jointSpaceInertia(mdl,Xup,Is_subt,H);

		copyMat(nv,nv,H,H_prev);
	} else {
		copyMat(nv,nv,H_prev,H);
	}

	//Cholesky decomposition of H
	if (!do_reuse_H || !H_init) {
		chol(nv,H,HL);
		copyMat(nv,nv,HL,HL_prev);
	} else {
		copyMat(nv,nv,HL_prev,HL);
	}
	H_init = true;

	//joint space bias force, C

	if (!do_approx_C) {
		jointSpaceBiasForce(mdl,Xup,qvel,C);

	} else {
		//approximate, all inertia at body frame

		setVec(nv,0.0,C);

		//spatial accel of gravity in body coords
		Vec6b a_grav;
		setcVec3(0.0, a_grav);
		mulcVec3(Xup[0]+COL2, mdl.grav, a_grav+COL1); //world z axis in body coords * grav acc

		Vec6b tmp;

		multMatVec6b(Is_subt[0],qvel,tmp);
		crossVec6bForce(qvel,tmp,C);
		multMatVec6b(Is_subt[0],a_grav,tmp);
		addVec6b(C,tmp,C);
	}
}

//u not const because u.err is output
void forwardDyn(const WmrModel& mdl, const Real state0[], const Real qvel0[], ControllerIO& u, 
	const HomogeneousTransform HT_parent[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, const Real dt, //inputs
	Real qacc[]) { //outputs
	//TODO, additional outputs class

	//for allocation
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//compute joint space inertia and bias force
	Real H[MAXNV*MAXNV];
	Real HL[MAXNV*MAXNV];
	Real C[MAXNV];

	HandC(mdl, HT_parent, qvel0, H, HL, C);

	if (mdl.use_constraints) {
		ConstraintJacobian A;

		constraintJacobians(mdl, state0, qvel0, HT_world, contacts, A );

		if (mdl.wheelGroundContactModel == 0) {
			forwardDynErpCfm(mdl, state0, qvel0, u.cmd, contacts, dt, HL, C, A, //input
				qacc); //output
		} else {
			forwardDynForceBalance(mdl, state0, qvel0, u, contacts, dt, HL, C, A, //input
				qacc); //output
		}
	} else {
		//TODO, forwardDynUnc()
	}
	
	return;

}

void constraintJacobians(const WmrModel& mdl, const Real state0[], const Real qvel0[], const HomogeneousTransform HT_world[], const ContactGeom* contacts, //input
	ConstraintJacobian& A ) { //output

	const int MAXNJ = WmrModel::MAXNF-1;

	//PRE-PROCESS

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nj = nf-1;
	const int nv = NUMQVEL(nf);

	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	const int na = mdl.get_na();

	const int* actframeinds = mdl.get_actframeinds();

	//END PRE-PROCESS

	//contact constraints
	int ncc = 0;
	if (nw > 0) {
		const WheelContactGeom* wcontacts;
		wcontacts = static_cast<const WheelContactGeom*>(contacts);
		ncc = wheelJacobians(mdl, HT_world, wcontacts, A.contact);
	} else if (nt > 0) {
		const TrackContactGeom* tcontacts;
		tcontacts = static_cast<const TrackContactGeom*>(contacts);
		ncc = trackJacobians(mdl, HT_world, tcontacts, A.contact);
	}
	

	//actuator constraints
	const int nac = na;
	if (nac > 0) {
		setMat(nac,nv,0.0,A.act);
		for (int ai=0; ai < nac; ai++) {
			int vi=TOQVELI(actframeinds[ai]);
			A.act[S2I(ai,vi,nac)] = 1.0;
		}
	}
	
	//(holonomic) joint constraints
	const int njc = mdl.get_njc();
	if (njc > 0) {

		Real c[WmrModel::MAXNJC]; //dummy
		Real Jc[WmrModel::MAXNJC*MAXNJ]; //njc x nj

		mdl.holonomicJointConstraints(mdl, state0+SI_JD, qvel0+VI_JR, //inputs
			c, Jc, 0, 0, 0); //outputs

		//set Ajoint
		setMatBlock(njc,0,0,njc,VI_JR,0.0, A.joint);
		copyMatBlock(njc,0,0,njc,nj,Jc, njc,0,VI_JR,A.joint);
	}

	//copy to A.all
	const int nc = ncc + nac + njc;
	const int contact_i0 = 0;
	const int act_i0 = ncc;
	const int joint_i0 = ncc + nac;

	if (ncc > 0)
		copyMatBlock(ncc,0,0,ncc,nv,A.contact, nc,contact_i0,0,A.all);
	if (nac > 0)
		copyMatBlock(nac,0,0,nac,nv,A.act, nc,act_i0,0,A.all);
	if (njc > 0)
		copyMatBlock(njc,0,0,njc,nv,A.joint, nc,joint_i0,0,A.all);
}

int parseWheelContacts(const WheelContactGeom* wcontacts, const int nw, bool incontact[], Real dz0[]) {
	int np = 0;
	np = nw;
	for (int wno=0; wno < nw; wno++) {
		incontact[wno] = wcontacts[wno].incontact;
		dz0[wno] = wcontacts[wno].dz;
	}
	return np;
}

int parseTrackContacts(const TrackContactGeom* tcontacts, const int nt,	bool incontact[], Real dz0[], int whichtrack[]) {
	int np = 0;
	for (int tno=0; tno < nt; tno++) {
		int np_ = tcontacts[tno].get_np(); //number of points for this track only
		copyVec(np_, tcontacts[tno].incontact, incontact+np);
		setVec(np_, tno, whichtrack+np);
		copyVec(np_, tcontacts[tno].dz, dz0+np);
		np += np_;
	}
	return np;
}


void forwardDynErpCfm(const WmrModel& mdl, const Real state0[], const Real qvel0[], const Real u_cmd[], const ContactGeom* contacts, const Real dt, 
	const Real HL[], const Real C[], const ConstraintJacobian& A, //input
	Real qacc[]) { //output
	//use erp, cfm like Open Dynamics Engine

	//for allocation
	const int MAXNJ = WmrModel::MAXNF-1;
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	const int MAXNP = ConstraintJacobian::MAXNP;
	const int MAXNC = ConstraintJacobian::MAXNC;

	// PRE-PROCESS
	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);
	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	const int na = mdl.get_na();

	const int* actframeinds = mdl.get_actframeinds();

	//contacts
	int np; //total number of points
	bool incontact[MAXNP];
	Real dz0[MAXNP]; //contact height errors
	int whichtrack[MAXNP];
	
	const WheelContactGeom* wcontacts;
	const TrackContactGeom* tcontacts;

	if (nw > 0) {
		wcontacts = static_cast<const WheelContactGeom*>(contacts);
		np = parseWheelContacts(wcontacts, nw, incontact, dz0);
	} else if (nt > 0) {
		tcontacts = static_cast<const TrackContactGeom*>(contacts);
		np = parseTrackContacts(tcontacts, nt, incontact, dz0, whichtrack);
	}

	assert(np <= MAXNP); //DEBUGGING

	int incontactinds[MAXNP];
	int npic; //number of points in contact
	npic = logicalFind(np, incontact, incontactinds);

	//count constraints
	const int ncc = 3*npic;
	const int nac = na;
	const int njc = mdl.get_njc();
	const int nc = ncc + nac + njc;

	//for constraint indexing
	const int contact_i0 = 0;
	const int act_i0 = ncc;
	const int joint_i0 = ncc + nac;

	// END PRE-PROCESS


	Real b[MAXNC];
	Real cfm_diag[MAXNC];
		
	if (ncc > 0) {

		Real vc0_incontact[3*MAXNP]; //contact point velocities for points in contact, concatenated into single vector
		multMatVec(3*npic,nv,A.contact,qvel0,1.0,vc0_incontact);

		for (int i=0; i<npic; i++) {
			int pno = incontactinds[i]; //point number in list of all
			int wtno; //wheel or track number
			if (nw > 0)
				wtno = pno;
			else if (nt > 0)
				wtno = whichtrack[pno];

			b[contact_i0+(i*3)+0] = - vc0_incontact[(i*3)+0];
			b[contact_i0+(i*3)+1] = - vc0_incontact[(i*3)+1];
			b[contact_i0+(i*3)+2] = -(mdl.erp_z[wtno] * dz0[pno])/dt - vc0_incontact[(i*3)+2];

			cfm_diag[contact_i0+(i*3)+0] = mdl.fds_x[wtno];
			cfm_diag[contact_i0+(i*3)+1] = mdl.fds_y[wtno];
			cfm_diag[contact_i0+(i*3)+2] = mdl.cfm_z[wtno];

		}
	}
		
	if (nac > 0) {

		for (int i=0; i<nac; i++) {
			int vi = TOQVELI(actframeinds[i]);

			b[act_i0+i] = u_cmd[i] - qvel0[vi];
			cfm_diag[act_i0+i] = 0;
		}
	}

	if (njc > 0) {

		Real c[WmrModel::MAXNJC];
		Real Jc[WmrModel::MAXNJC*MAXNJ]; //dummy

		mdl.holonomicJointConstraints(mdl, state0+SI_JD, qvel0+VI_JR, //inputs
			c, Jc, 0, 0, 0); //outputs

		Real cd0[WmrModel::MAXNJC]; //d/dt c based on qvel0
		multMatVec(njc,nv,A.joint,qvel0,1.0,cd0);

		for (int i=0; i<njc; i++) {
			b[joint_i0+i] = -mdl.erp_j[i]*c[i]/dt - cd0[i];
			cfm_diag[joint_i0+i] = mdl.cfm_j[i];
		}
	}
	mulcVec(nc,1/dt,b);
	
	//TODO, duplication

	//applied force
	Real tau[MAXNV];
	setVec(nv,0.0,tau);

	Real tau_minus_C[MAXNV]; //tau - C
	//TODO, make this a single operation?
	copyVec(nv,tau,tau_minus_C);
	addmVec(nv,C,-1.0,tau_minus_C);

	Real invH_tau_minus_C[MAXNV]; //inv(H)*(tau-C)
	cholSolve(nv, HL, tau_minus_C, invH_tau_minus_C);

	//end duplication

	//Hl*lambda = fl
		
	Real fl[MAXNC];
	//fl = b - A*inv(H)*(tau-C)

	multMatVec(nc, nv, A.all, invH_tau_minus_C, -1.0, fl);
	addmVec(nc, b, 1.0, fl);

	//Hl = A*inv(H)*A^T + CFM
	Real Hl[MAXNC*MAXNC];
	Real AT[MAXNV*MAXNC]; //A^T
	copyTMat(nc,nv,A.all,AT); //TODO, eliminate this?
	Real invH_AT[MAXNV*MAXNC]; //inv(H)*A^T
	cholSolveMat(nv,HL,nc,AT,invH_AT);
	multMatMat(nc,nv,A.all,nc,invH_AT,1.0,Hl);

	//add cfm
	for (int i=0; i<nc; i++)
		Hl[S2I(i,i,nc)] += cfm_diag[i]/dt;

	//DEBUGGING
	//std::cout << "Hl=\n"; printMatReal(nc,nc,Hl,-1,-1);
	//std::cout << "fl=\n"; printMatReal(nc,1,fl,-1,-1);

	//Cholesky decomposition of Hl
	Real HlL[MAXNC*MAXNC]; //lower triangular
	chol(nc,Hl,HlL);


	Real lambda[MAXNC]; //constraint forces
	cholSolve(nc, HlL, fl, lambda);

	//compute acceleration
	//qacc = inv(H)*(tau-C + A^T*lambda)
	//qacc = inv(H)*(tau-C) + (inv(H)*A^T)*lambda
	multMatVec(nv,nc,invH_AT,lambda,1.0,qacc);
	addmVec(nv,invH_tau_minus_C,1.0,qacc);

	return;
}


void forwardDynForceBalance(const WmrModel& mdl, const Real state0[], const Real qvel0[], ControllerIO& u, const ContactGeom* contacts, const Real dt, 
	const Real HL[], const Real C[], const ConstraintJacobian& A, //input
	Real qacc[]) { //output

	//OPTIONS
	const Real qr_tol = 1e-3;
	const int max_iter = 20;
	const Real cost_tol = 1;
	const Real dcost_tol = cost_tol/10;


	//for allocation
	const int MAXNJ = WmrModel::MAXNF-1;
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	const int MAXNP = ConstraintJacobian::MAXNP;
	const int MAXNC = ConstraintJacobian::MAXNC;

	// PRE-PROCESS
	const int nf = mdl.get_nf();
	const int nj = nf - 1;
	const int nv = NUMQVEL(nf);
	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	const int na = mdl.get_na();

	const int* wheelframeinds = mdl.get_wheelframeinds();
	const int* sprocketframeinds = mdl.get_sprocketframeinds();
	const int* actframeinds = mdl.get_actframeinds();
	const Frame* frames = mdl.get_frames();

	//contacts
	int np; //total number of points
	bool incontact[MAXNP];
	Real dz0[MAXNP]; //contact height errors
	int whichtrack[MAXNP];

	const WheelContactGeom* wcontacts;
	const TrackContactGeom* tcontacts;

	if (nw > 0) {
		wcontacts = static_cast<const WheelContactGeom*>(contacts);
		np = parseWheelContacts(wcontacts, nw, incontact, dz0);
	} else if (nt > 0) {
		tcontacts = static_cast<const TrackContactGeom*>(contacts);
		np = parseTrackContacts(tcontacts, nt, incontact, dz0, whichtrack);
	}

	assert(np <= MAXNP); //DEBUGGING

	int incontactinds[MAXNP];
	int npic; //number of points in contact
	npic = logicalFind(np, incontact, incontactinds);

	//count constraints
	const int ncc = 3*npic;
	const int nac = na;
	const int njc = mdl.get_njc();
	const int nc = ncc + nac + njc;

	//for constraint indexing
	const int contact_i0 = 0;
	const int act_i0 = ncc;
	const int joint_i0 = ncc + nac;

	// END PRE-PROCESS

	//initialize inputs
	Real xall[MAXNC]; //size nc
	Real x[MAXNC]; //size ninpt


	//contact constraints
	Real vc0_incontact[3*MAXNP]; //contact point velocities for points in contact, concatenated into single vector
	if (ncc > 0) {

		multMatVec(3*npic,nv,A.contact,qvel0,1.0,vc0_incontact);

		//init to contact pt velocity at previous time step
		copyVec(ncc,vc0_incontact,xall+contact_i0);

		//init vz to zero, or never reaches steady state
		for (int i=0; i<npic; i++) 
			xall[contact_i0+(3*i)+2] = 0;
	}

	Real u_act0[WmrModel::MAXNA]; //na x 1, from qvel0
	if (nac > 0) {

		for (int ai = 0; ai < na; ai++) {
			int vi = TOQVELI(actframeinds[ai]);
			u_act0[ai] = qvel0[vi];
		}

		copyVec(nac,u_act0, xall+act_i0); //init to joint velocity at previous time step
	}

	
	Real cd0[WmrModel::MAXNJC]; //d/dt c based on qvel0
	if (njc > 0) {

		multMatVec(njc,nv,A.joint,qvel0,1.0,cd0);

		copyVec(njc,cd0,xall+joint_i0);
	}

	//TODO, move to subfunction?
	//get maximal linearly independent subset of constraints (rows of A)
	//set includes all actuator constraints, some wheel and joint constraints
	
	struct ConstraintBool {
		bool ind[MAXNC]; //independent
		bool inpt[MAXNC]; //input
	};
	ConstraintBool cis;


	setVec(nc, true, cis.ind); //initialize
	
	int nind; //number of constraints in independent subset
	int ncc_ind;
	int njc_ind;

	
	if (nc-nac > 0) {

		//is NOT actuator constraint
		bool cis_nact[MAXNC];
		setVec(nc,true,cis_nact);
		setVec(nac,false,cis_nact+act_i0);

		//is actuated DOF in joint space vel
		bool vis_act[MAXNV];
		setVec(VI_JR,false,vis_act);
		for (int fi = 1; fi < nf; fi++)
			vis_act[TOQVELI(fi)] = frames[fi].isactuated;

		Real A_[MAXNC*MAXNV]; //submatrix of A, (nc-nac) x (nv-na)
		int col = 0;
		for (int vi=0; vi < nv; vi++) {
			if (!vis_act[vi]) {
				logicalIndexIn(nc,cis_nact,A.all+S2I(0,vi,nc),A_+S2I(0,col,nc-nac));
				col++;
			}
		}

		//DEBUGGING
		//std::cout << "A=\n"; printMatReal(nc,nv,A,-1,-1);
		//std::cout << "A_=\n"; printMatReal(nc-nac,nv-na,A_,-1,-1);

		bool is_ind_[MAXNC];
		subset(nc-nac, nv-na, A_, qr_tol, is_ind_);

		logicalIndexOut(nc,cis_nact,is_ind_,cis.ind);
		ncc_ind = logicalCount(ncc,cis.ind + contact_i0);
		njc_ind = logicalCount(njc,cis.ind + joint_i0);
		nind = ncc_ind + nac + njc_ind;

	} else {
		nind = nc;
		ncc_ind = ncc;
		njc_ind = njc;
	}

	//for indexing the independent subset of constraints
	int contact_i0_ind, act_i0_ind, joint_i0_ind;
	contact_i0_ind = 0;
	act_i0_ind = ncc_ind;
	joint_i0_ind = ncc_ind + nac;


	//eliminate dependent rows from A
	const int MAXNIND = MAXNV; //max number of independent constraints, = max number of inputs

	Real Aind[MAXNIND*MAXNV]; 
	int row=0; //row in Aind
	for (int ci=0; ci<nc; ci++) {
		if (cis.ind[ci]) {
			copyMatRow(nc,nv,ci,A.all, nind,row,Aind);
			row++;
		}
	}

	//DEBUGGING
	//std::cout << "Aind=\n"; printMatReal(nind,nv,Aind,-1,-1);

	//independent wheel, joint constraints
	Real Acontact_ind[MAXNIND*MAXNV]; 
	Real Ajoint_ind[WmrModel::MAXNJC*MAXNV];

	copyMatBlock(nind,contact_i0_ind,0,ncc_ind,nv,Aind, ncc_ind,0,0,Acontact_ind);
	copyMatBlock(nind,joint_i0_ind,0,njc_ind,nv,Aind, njc_ind,0,0,Ajoint_ind);


	//is input
	copyVec(nc, cis.ind, cis.inpt);
	if (nac > 0 && mdl.actuatorModel == 0) //ideal actuators
		setVec(nac, false, cis.inpt + act_i0);
	int ninpt = logicalIndexIn(nc,cis.inpt,xall,x);


	//PRECOMPUTE VARS FOR COST FUNCTION

	//TODO, duplication
	//applied force
	Real tau[MAXNV];
	setVec(nv,0.0,tau);

	Real tau_minus_C[MAXNV]; //tau - C
	//TODO, make this a single operation?
	copyVec(nv,tau,tau_minus_C);
	addmVec(nv,C,-1.0,tau_minus_C);

	Real invH_tau_minus_C[MAXNV]; //inv(H)*(tau-C)
	cholSolve(nv, HL, tau_minus_C, invH_tau_minus_C);

	//end duplication

	//Hl*lambda = fl

	Real fl_c[MAXNV]; //constant part of fl
	//fl_c = Aind*inv(H)*(tau-C)
	multMatVec(nind, nv, Aind, invH_tau_minus_C, -1.0, fl_c);

	//Hl = Aind*inv(H)*Aind^T
	Real Hl[MAXNIND*MAXNIND];
	Real AindT[MAXNV*MAXNIND]; //Aind^T
	copyTMat(nind,nv,Aind,AindT); //TODO, eliminate this?
	Real invH_AindT[MAXNV*MAXNIND]; //inv(H)*Aind^T
	cholSolveMat(nv,HL,nind,AindT,invH_AindT);
	multMatMat(nind,nv,Aind,nind,invH_AindT,1.0,Hl);


	//Cholesky decomposition of Hl
	Real HlL[MAXNIND*MAXNIND]; //lower triangular
	chol(nind,Hl,HlL);
	
	Real cost;

	//additional outputs of cost function
	Real mf[MAXNC]; //model force

	struct dModelForce {
		Real contact_dwgcin[MAXNP][3*5];
		Real act_du[WmrModel::MAXNA*WmrModel::MAXNA];
		Real joint_djd[WmrModel::MAXNJC*MAXNJ];
		Real joint_djr[WmrModel::MAXNJC*MAXNJ];
	};
	dModelForce dmf;

	Real err[MAXNV]; //force-balance error

	//calculate cost
	
	//lambda closure, requires C++11
	auto forwardDynCost = [&] ( const Real x_[] ) mutable -> Real {
		
		//the following must be mutable: qacc, mf, dmf, err

		//why must redefine constants?
		const int MAXNJ = WmrModel::MAXNF-1;
		const int MAXNV = NUMQVEL(WmrModel::MAXNF);
		const int MAXNP = ConstraintJacobian::MAXNP;
		const int MAXNC = ConstraintJacobian::MAXNC;
		const int MAXNIND = MAXNV;


		//x_ is size ninpt
		Real xb[MAXNC]; //buffered, size nc
		logicalIndexOut(nc,cis.inpt,x_,xb);

		Real b[MAXNC]; //size nc

		if (ncc > 0) {
			copyVec(ncc, xb+contact_i0, b+contact_i0);
			addmVec(ncc, vc0_incontact, -1.0, b+contact_i0);
		}
		if (nac > 0) {
			if (mdl.actuatorModel != 0)	
				copyVec(nac, xb+act_i0, b+act_i0);
			else 
				copyVec(nac, u.cmd, b+act_i0); //ideal actuators
			addmVec(nac, u_act0, -1.0, b+act_i0);
		}
		if (njc > 0) {
			copyVec(njc, xb+joint_i0, b+joint_i0);
			addmVec(njc, cd0, -1.0, b+joint_i0);
		}

		Real b_[MAXNIND]; //size nind
		logicalIndexIn(nc, cis.ind, b, b_);

		mulcVec(nind, 1/dt, b_);

		//fl = b + fl_c
		Real fl[MAXNIND];
		copyVec(nind, b_, fl);
		addmVec(nind, fl_c, 1.0, fl);


		Real lambda[MAXNIND]; //constraint forces
		cholSolve(nind, HlL, fl, lambda);
		
		//DEBUGGING
		//std::cout << "lambda=\n"; printMatReal(nind,1,lambda,-1,-1);

		//compute acceleration
		//qacc = inv(H)*(tau-C + Aind^T*lambda)
		//qacc = inv(H)*(tau-C) + (inv(H)*Aind^T)*lambda
		multMatVec(nv,nind,invH_AindT,lambda,1.0,qacc);
		addmVec(nv,invH_tau_minus_C,1.0,qacc);

		Real qvel[MAXNV]; //new qvel, at time i+1
		copyVec(nv,qvel0,qvel);
		addmVec(nv,qacc,dt,qvel);

		//force-balance error
		
		if (ncc > 0) {
			//inputs for wheel-ground contact model, time i+1
			//only for points in contact!
			Vec3 vc[MAXNP]; //velocity of contact pts
			Real Rw[MAXNP]; //wheel radius * wheel angular rate
			Real dz[MAXNP]; //contact height error

			Vec3 fc[MAXNP]; //contact force (output)

			Real tmp[3*MAXNP]; //temporary
			multMatVec(ncc,nv,A.contact,qvel,1.0,tmp);
			copyArrayToVec3(tmp,npic,vc);
			
			if (nw > 0) {
				for (int i=0; i<npic; i++) { //loop over points in contact
					int wno = incontactinds[i];
					int fi = wheelframeinds[wno];
					Rw[i] = frames[fi].rad * qvel[TOQVELI(fi)];
					dz[i] = wcontacts[wno].dz + vc[i][2]*dt; // + vz*dt, symplectic Euler

					mdl.wheelGroundContactModel(wno, mdl.wgc_p, vc[i], Rw[i], dz[i], //inputs
						fc[i], (Real*) dmf.contact_dwgcin+(i*3*5)); //outputs
				}


			} else if (nt > 0) {
				for (int i=0; i<npic; i++) {
					int pno = incontactinds[i];
					int tno = whichtrack[pno];
					int fi = sprocketframeinds[tno];

					Rw[i] = frames[fi].rad * qvel[TOQVELI(fi)];
					dz[i] = dz0[pno] + vc[i][2]*dt; // + vz*dt, symplectic Euler

					mdl.wheelGroundContactModel(tno, mdl.wgc_p, vc[i], Rw[i], dz[i], //inputs
						fc[i], (Real*) dmf.contact_dwgcin+(i*3*5)); //outputs
				}
			}

			copyVec3ToArray(npic, fc, mf+contact_i0);
		}
		if (nac > 0 && mdl.actuatorModel != 0) {
			Real u_act[WmrModel::MAXNA];
			for (int ai=0; ai < na; ai++)
				u_act[ai] = qvel[TOQVELI(actframeinds[ai])];

			mdl.actuatorModel( mdl.act_p, u.cmd, u_act, u.interr, mf+act_i0, u.err, dmf.act_du);
		}
		if (njc > 0) {
			Real jd[MAXNJ]; //joint displacement at next time step
			copyVec(nj, state0+SI_JD, jd);
			addmVec(nj, qvel+VI_JR, dt, jd); //symplectic Euler

			//dummy
			Real c[WmrModel::MAXNJC];
			Real Jc[WmrModel::MAXNJC * MAXNJ];

			mdl.holonomicJointConstraints(mdl, jd, qvel+VI_JR, c, Jc, mf+joint_i0, dmf.joint_djd, dmf.joint_djr);

		}

		Real tau_lambda[MAXNV];
		Real tau_mf[MAXNV];

		multMatTVec(nind,nv,Aind,lambda,1.0,tau_lambda);
		multMatTVec(nc,nv,A.all,mf,1.0,tau_mf);

		for (int vi=0; vi < nv; vi++) 
			err[vi] = tau_lambda[vi] - tau_mf[vi];

		if (nac > 0 && mdl.actuatorModel == 0) { //ideal actuators
			for (int ai=0; ai<nac; ai++)
				err[TOQVELI(actframeinds[ai])] = 0.0;
		}

		Real cost = dotVec(nv,err,err);

		//DEBUGGING
		//std::cout << "tau_lambda=\n"; printMatReal(nv,1,tau_lambda,-1,-1);
		//std::cout << "tau_mf=\n"; printMatReal(nv,1,tau_mf,-1,-1);
		//std::cout << "err=\n"; printMatReal(nv,1,err,-1,-1);

		return cost;

	};
	cost = forwardDynCost(x);


	//precomputed vars
	Real dwgcin_dx[MAXNP][5*MAXNIND]; //nwic x 5 x ninpt, Jacobian of [vx vy vz Rw dz] wrt x
	Real du_dx[WmrModel::MAXNA*MAXNIND]; //na x ninpt
	Real djd_dx[MAXNJ*MAXNIND]; //nj x ninpt
	Real djr_dx[MAXNJ*MAXNIND]; //nj x ninpt
	Real dtaulambda_dx[MAXNV*MAXNIND]; //nv x ninpt

	//additional outputs of gradient function
	Real derr_dx[MAXNV*MAXNIND];
	Real grad[MAXNIND]; // 1 x ninpt, dcost/dx

	int iter=0;

	//std::cout << "iter= " << iter << ", cost= " << cost << std::endl; //DEBUGGING

	if (cost > cost_tol) { //only do this if necessary

		//PRECOMPUTE VARS FOR GRADIENT CALCULATION
		auto forwardDynPrecompute = [&] ( void ) mutable {
			//the following must be mutable: dwgcin_dx, du_dx, djd_dx, djr_dx, dtaulambda_dx

			const int MAXNV = NUMQVEL(WmrModel::MAXNF);
			const int MAXNIND = MAXNV;

			Real dlambda_dx[MAXNIND*MAXNIND]; //nind x ninpt

			Real inv_Hl[MAXNIND*MAXNIND]; //TODO, avoid computing this?
			Real I_nind[MAXNIND*MAXNIND]; //identity matrix
			setMat(nind,nind,0.0,I_nind);
			for (int i=0; i<nind; i++) 
				I_nind[S2I(i,i,nind)]=1.0;
			
			cholSolveMat(nind,HlL,nind,I_nind,inv_Hl);

			//dlambda_dx, a subset of columns of inv_Hl
			{
				int indno=0; //column index in Hl
				int inptno=0; //column index in dlambda_dx;
				for (int i=0; i<nc; i++) { //loop over constraints
					if (cis.ind[i]) {
						if (cis.inpt[i]) {
							copyMatCol(nind,indno,inv_Hl, inptno,dlambda_dx);
							inptno++;
						}
						indno++;
					}
				}
			}


			mulcMat(nind,ninpt,1/dt,dlambda_dx);

			//dqvel/dx = inv(H)*Aind^T * dlambda/dx (*dt)
			Real dqvel_dx[MAXNV*MAXNIND]; //nv x ninpt
			multMatMat(nv,nind,invH_AindT,ninpt,dlambda_dx,dt,dqvel_dx);
		

			//std::cout << "dlambda_dx=\n"; printMatReal(nv,ninpt,dlambda_dx,-1,-1);
			//std::cout << "dqvel_dx=\n"; printMatReal(nv,ninpt,dqvel_dx,-1,-1);
		

			if (ncc > 0) {
				Real Acontact_[3*MAXNV]; // submatrix for a single point
				Real dvc_dx[3*MAXNIND]; // 3 x ninpt, temporary

				for (int i=0; i<npic; i++) {
					int fi; //frame index
					if (nw > 0) {
						int wno = incontactinds[i];
						fi = wheelframeinds[wno];
					} else if (nt > 0) {
						int tno = whichtrack[incontactinds[i]];
						fi = sprocketframeinds[tno];
					}

					//d [vx vy vz] /dx
					//TODO, eliminate copy?
					copyMatBlock(ncc,contact_i0+3*i,0,3,nv,A.contact, 3,0,0,Acontact_);
					multMatMat(3,nv,Acontact_,ninpt,dqvel_dx,1.0,dvc_dx);
					copyMatBlock(3,0,0,3,ninpt,dvc_dx, 5,0,0,dwgcin_dx[i]);
					//d Rw /dx
					copyMatRow(nv,ninpt,TOQVELI(fi),dqvel_dx, 5,3,dwgcin_dx[i]);
					mulcMatRow(5,ninpt,3,frames[fi].rad,dwgcin_dx[i]);
					//d dz /dx = dvz/dx*dt
					copyMatRow(5,ninpt,2,dwgcin_dx[i], 5,4,dwgcin_dx[i]);
					mulcMatRow(5,ninpt,4,dt,dwgcin_dx[i]);

				}

			}
			if (nac > 0 && mdl.actuatorModel != 0) {
				//du/dx
				for (int ai=0; ai<nac; ai++) {
					int vi=TOQVELI(actframeinds[ai]);
					copyMatRow(nv,ninpt,vi,dqvel_dx, nac,ai,du_dx);
				}
			}
			if (njc > 0) {
				// djr/dx and djd/dx (= djr/dx*dt)
				copyMatBlock(nv,VI_JR,0,nj,ninpt,dqvel_dx, nj,0,0,djr_dx);
				copyMat(nj,ninpt,djr_dx,djd_dx);
				mulcMat(nj,ninpt,dt,djd_dx);
			}
			multMatTMat(nind,nv,Aind,ninpt,dlambda_dx,1.0,dtaulambda_dx);
		};
		forwardDynPrecompute();

		//optimize
		Real cost_prev = REALMAX;
		Real Hess[MAXNIND*MAXNIND]; //ninpt x ninpt
		Real HessL[MAXNIND*MAXNIND]; //lower triangular
		Real p[MAXNIND];

		
		//calculate gradient

		//lambda closure, requires C++11
		auto forwardDynGradient = [&] ( Real grad_[] ) mutable {

			//the following must be mutable: derr_dx

			//why must redefine constants?
			const int MAXNV = TOQVELI(WmrModel::MAXNF);
			const int MAXNP = ConstraintJacobian::MAXNP;
			const int MAXNC = ConstraintJacobian::MAXNC;
			const int MAXNIND = MAXNV;

			Real dmf_dx[MAXNC*MAXNIND]; //nc x ninpt

			if (ncc > 0) {
				Real dmfcontact_dx[3*MAXNIND]; // 3 x ninpt, for a single point
				
				for (int i=0; i<npic; i++) { //loop over points in contact
					multMatMat(3,5,(Real*)dmf.contact_dwgcin+(i*3*5), ninpt,(Real*)dwgcin_dx+(i*5*MAXNIND), 1.0, dmfcontact_dx);
					copyMatBlock(3,0,0,3,ninpt,dmfcontact_dx, nc,contact_i0+(i*3),0,dmf_dx);
				}
			}
			if (nac > 0 && mdl.actuatorModel != 0) {
				Real dmfact_dx[WmrModel::MAXNA*MAXNIND]; //na x ninpt
				multMatMat(nac,nac,dmf.act_du,ninpt,du_dx,1.0,dmfact_dx);
				copyMatBlock(na,0,0,na,ninpt,dmfact_dx, nc,act_i0,0,dmf_dx);
			}
			if (njc > 0) {
				Real dmfjoint_dx[WmrModel::MAXNJC*MAXNIND]; //nv x ninpt

				multMatMat(njc,nj,dmf.joint_djd,ninpt,djd_dx,1.0,dmfjoint_dx);
				copyMatBlock(njc,0,0,njc,ninpt,dmfjoint_dx, nc,joint_i0,0,dmf_dx); //dmf_joint_djd * djd_dx

				multMatMat(njc,nj,dmf.joint_djr,ninpt,djr_dx,1.0,dmfjoint_dx); //dmf_joint_djr * djr_dx
				addmMatBlock(njc,0,0,njc,ninpt,dmfjoint_dx, nc,joint_i0,0,1.0,dmf_dx);
			}

			Real dtaumf_dx[MAXNV*MAXNIND]; //nv x ninpt
			multMatTMat(nc,nv,A.all,ninpt,dmf_dx,1.0,dtaumf_dx); //A*dtaumf/dx

			for (int i=0; i<nv*ninpt; i++) 
				derr_dx[i] = dtaulambda_dx[i] - dtaumf_dx[i];


			if (nac > 0 && mdl.actuatorModel == 0) { //ideal actuators
				for (int ai=0; ai<nac; ai++) {
					int fi = actframeinds[ai];
					setMatRow(nv,ninpt,TOQVELI(fi),0.0,derr_dx);
				}
			}

			multMatTMat(nv,1,err,ninpt,derr_dx,2.0,grad_); //grad = 2*err^T*(derr/dx)

			//DEBUGGING
			//std::cout << "dtaulambda_dx=\n"; printMatReal(nv,ninpt,dtaulambda_dx,-1,-1);
			//std::cout << "dtaumf_dx=\n"; printMatReal(nv,ninpt,dtaumf_dx,-1,-1);
			//std::cout << "derr_dx=\n"; printMatReal(nv,ninpt,derr_dx,-1,-1);
			//std::cout << "grad_=\n"; printMatReal(1,ninpt,grad_,-1,-1);

			return;

		};
		forwardDynGradient(grad);
		

		for (iter=1; iter<max_iter; iter++) {

			//compute Hessian, H = 2*(drre_dx^T * derr_dx)
			multMatTMat(nv,ninpt,derr_dx,ninpt,derr_dx,2.0,Hess);

			chol(ninpt,Hess,HessL);
			cholSolve(ninpt,HessL,grad,p);

			mulcVec(ninpt,-1.0,p); //negate

			Real gradp = dotVec(ninpt,grad,p);
			if (fabs(gradp) < dcost_tol) {
				//gradient is zero at minimum, stop if close enough
				break;
			}

			Real alpha;

			alpha = linesearch( ninpt, p, 1.0, forwardDynCost, forwardDynGradient, x, cost, grad);

			//std::cout << "Hess = \n"; printMatReal(ninpt,ninpt,Hess,-1,-1);
			//std::cout << "p = \n"; printMatReal(ninpt,1,p,-1,-1);

			//DEBUGGING, force-balance stabilization, take a single Newton's method step
			//if (0)
			//	break;

			//std::cout << "iter= " << iter << ", cost= " << cost << std::endl; //DEBUGGING

			if (cost < cost_tol) 
				break;
			if (fabs(cost - cost_prev) < dcost_tol) 
				break;
			cost_prev = cost;

		}
		
		
	}	
	
	return;
}


void odeDyn(const Real time, const Real y[], const WmrModel& mdl, const SurfaceVector& surfaces, ContactGeom* contacts, const Real dt, //inputs
	Real ydot[], HomogeneousTransform HT_parent[] ) { //outputs

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int na = mdl.get_na();

	const int ns = NUMSTATE(nf);
	const int nv = NUMQVEL(nf);
	const Frame* frames = mdl.get_frames();
	
	ControllerIO u;

	//parse y
	Real state[MAXNS];
	Real qvel[MAXNV];

	copyVec(ns,y,state);
	copyVec(nv,y+ns,qvel);
	copyVec(na,y+ns+nv,u.interr);

	//get control inputs
	mdl.controller(mdl, time, y, u.cmd, 0);

	//convert state to Homogeneous Transforms
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,y,HT_parent,HT_world);

	//update contact geometry
	updateModelContactGeom(mdl, surfaces, HT_world, 0, contacts);

	//compute joint space acceleration
	Real qacc[MAXNV];

	forwardDyn(mdl, state, qvel, u, HT_parent, HT_world, contacts, dt, qacc);

	if (dt==dt) { //not NaN
		//Semi-implicit or symplectic Euler, conserves energy
		//http://en.wikipedia.org/wiki/Semi-implicit_Euler_method
		//v(i+1) = v(i) + g(x(i),v(i))*dt
		//x(i+1) = x(i) + f(x(i),v(i+1))*dt
		//standard Euler uses v(i) to compute x(i+1) instead of v(i+1)

		//given dt, forwardDyn.m assumes symplectic Euler
		//integrate qacc to update qvel
		addmVec(nv,qacc,dt,qvel);
	}

	//convert qvel to time derivative of state
	Real statedot[MAXNS];
	qvelToQdot(nf,qvel,state+SI_ORIENT,HT_world[0],statedot);

	//copy to ydot
	copyVec(ns,statedot,ydot);
	copyVec(nv,qacc,ydot+ns);
	copyVec(na,u.err,ydot+ns+nv);

}
