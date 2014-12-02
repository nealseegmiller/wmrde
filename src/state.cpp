#include<state.h>

//HT_parent[i] = transform from frame i coords to parent frame coords
//HT_world[i] = transform from frame i coords to world coords
void stateToHT(const WmrModel &mdl, const Real state[], HomogeneousTransform HT_parent[], HomogeneousTransform HT_world[]) {
	
	int nf = mdl.get_nf(); //number of frames
	const Frame* frames = mdl.get_frames();

	poseToHT(state+SI_ORIENT,state+SI_POS,HT_parent[0]); //body frame

	for (int fi=1; fi<nf; fi++) { //frame index
		Real jd = state[TOSTATEI(fi)];
		int dof_type = frames[fi].dof_type;
		const Real* HT_parent_jd0 = frames[fi].HT_parent_jd0;

		if (dof_type < 3) { //rotary
			Mat3 R;
			RotInd(dof_type, jd, R);
			multMatMat3(HT_parent_jd0, R, HT_parent[fi]); //compute rotation
			copyVec3(HT_parent_jd0+COL3, HT_parent[fi]+COL3); //copy translation
		} 
		else { //prismatic
			Vec3 t;
			copyRowToVec3(dof_type-3, HT_parent_jd0, t);
			mulcVec3(t, jd, t);
			addVec3(t, HT_parent_jd0+COL3, HT_parent[fi]+COL3); //compute translation
			copyMat3(HT_parent_jd0, HT_parent[fi]); //copy rotation
		}
	}

	if (HT_world != 0) { //if not null pointer
		copyHT(HT_parent[0],HT_world[0]);
		for (int fi=1; fi<nf; fi++) { //frame index
			int parent_ind = frames[fi].parent_ind;
			composeHT(HT_world[parent_ind], HT_parent[fi], HT_world[fi]);
		}
	}
}


//nf:		number of frames
//qvel:		joint space velocity
//orient:	orientation vector
//R_body_to_world:	Rotation matrix, redundant with orient, assumed to be consistent!
//qdot:		d/dt state
void qvelToQdot(const int nf, const Real qvel[], const Real orient[], const Mat3 R_body_to_world, Real qdot[]) {

	int nj = nf-1; //number of joints
	for (int ji=0; ji<nj; ji++) //joint index 
		qdot[ji+SI_JD] = qvel[ji+6]; //copy joint rates
	//TODO, unsafe if vectorized!
	multMatVec3(R_body_to_world, qvel+3, qdot+SI_POS); //rotate linear velocity to world coords
	velToOrientrate(orient, qvel+0, qdot+SI_ORIENT, 0); //convert angular velocity to orient rate
}

void qdotToQvel(const int nf, const Real qdot[], const Real orient[], const Mat3 R_body_to_world, Real qvel[]) {
	//this function is not necessary to simulate WMR

	int nj = nf-1;
	for (int ji=0; ji<nj; ji++) //joint index
		qvel[ji+6] = qdot[ji+SI_JD];
	multMatTVec3(R_body_to_world, qdot+SI_POS, qvel+3); //rotate linear velocity to body coords (transpose)
	orientrateToVel(orient, qdot+SI_ORIENT, qvel+0, 0); //convert orient rate to angular velocity
}

//compute the spatial velocity of each frame from the joint space velocity
//using the RNEA, compare to jointSpaceBiasForce
//Xup:		size nf array of Plucker transforms, Xup[i] transforms spatial *motion* vector from parent(i) to i coords
//qvel:		size nv, joint space velocity
//v:		size nf, spatial velocity of each frame
//TODO, check this
void qvelToSpatialVel(const WmrModel &mdl, Mat6b Xup[], const Real qvel[], Vec6b v[]) {

	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	const Frame* frames = mdl.get_frames();

	//body frame
	copyArrayToVec6b(qvel,v[0]);

	for (int fi=1; fi < nf; fi++) {
		int parent_fi = frames[fi].parent_ind;
		int dof_type = frames[fi].dof_type;

		Vec6b vJ; //spatial joint velocity
		setcVec6b(0.0,vJ);
		vJ[DofTypeToVec6bInd(dof_type)] = qvel[TOQVELI(fi)];

		//v = Xup * v(parent) + vJ
		multPluckerVec6b(Xup[fi], v[parent_fi], v[fi]);
		addVec6b(v[fi],vJ,v[fi]);
	}
}



