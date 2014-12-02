#include <WmrModel.h>
//WmrModel class implementation file

//constructor
Frame::Frame() {
	//do nothing
}

void Frame::setMass(const Real scalar_mass, const Vec3 center_of_mass, const Mat3 moment_of_inertia) {
	mass = scalar_mass;
	copyVec3(center_of_mass,cm);
	copyMat3(moment_of_inertia,I);
	toSpatialInertia(mass,cm,I,Is);
}

//static only in declaration, not definition
int Frame::convertDofTypeString(std::string dof_string) {
	int dof_int;
	if		( dof_string.compare("RX")==0 ) //revolute
		dof_int=0;
	else if ( dof_string.compare("RY")==0 ) 
		dof_int=1;
	else if ( dof_string.compare("RZ")==0 ) 
		dof_int=2;
	else if ( dof_string.compare("PX")==0 ) //prismatic
		dof_int=3;
	else if ( dof_string.compare("PY")==0 ) 
		dof_int=4;
	else if ( dof_string.compare("PZ")==0 ) 
		dof_int=5;
	else
		dof_int=-1; //invalid
	return dof_int;
}

//constructor
WmrModel::WmrModel() {
	nf=0;
	nw=0;
	na=0;
}

void WmrModel::initFrame(const int fi) {
	//initialize frame
	frames[fi].dof_type = -1;
	frames[fi].parent_ind = -1;
	setIdentityHT(frames[fi].HT_parent_jd0);
	
	frames[fi].iswheel = false;
	frames[fi].issprocket = false;
	frames[fi].isactuated = false;
	frames[fi].isfixed = false;

	//no mass
	Real scalar_mass = 0;
	Vec3 center_of_mass = {0,0,0};
	Mat3 moment_of_inertia;
	setcMat3(0.0, moment_of_inertia);

	frames[fi].setMass(scalar_mass, center_of_mass, moment_of_inertia);
}

//TODO, use logicalFind()?
void WmrModel::update_wheelframeinds() {
	nw = 0; //also updates nw
	for (int fi = 0; fi < nf; fi++) {
		if (frames[fi].iswheel) {
			wheelframeinds[nw] = fi;
			nw++;
		}
	}
}

void WmrModel::update_sprocketframeinds() {
	nt = 0; //also updates nt
	for (int fi = 0; fi < nf; fi++) {
		if (frames[fi].issprocket) {
			sprocketframeinds[nt] = fi;
			nt++;
		}
	}
}

void WmrModel::update_actframeinds() {
	na = 0; //also updates na
	for (int fi = 0; fi < nf; fi++) {
		if (frames[fi].isactuated) {
			actframeinds[na] = fi;
			na++;
		}
	}
}


void WmrModel::addBodyFrame(const std::string name) {

	assert(nf==0);

	//don't need this-> ?
	int i = 0;
	nf++;
	initFrame(i);
	frames[i].name = name;
	
	update_wheelframeinds(); 
	update_actframeinds();

	njc=0;
	grav=9.81; //m/s^2, earth
	//grav=1.62; //m/s^2, moon

	use_constraints = true;

	//init function pointers to zero
	wheelGroundContactModel = 0;
	actuatorModel = 0;
	holonomicJointConstraints = 0;
	controller = 0;

}


void WmrModel::addJointFrame(const std::string name, const std::string parent_name, const std::string dof_string, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0) {
	
	//error checking
	int mi; //match index
	assert(nf+1 <= MAXNF); //must not exceed allocated space

	mi = nameToInd(name);
	assert(mi<0); //name input must not match name of existing frame

	mi = nameToInd(parent_name);
	assert(mi >= 0); //parent_name must match name of existing frame
	assert(!frames[mi].iswheel); //parent can not be wheel frame

	int dof_int = Frame::convertDofTypeString(dof_string);
	assert(dof_int >= 0);
	
	//end error checking

	int i = nf;
	nf++;

	initFrame(i);
	frames[i].name = name;
	frames[i].parent_ind = mi;
	frames[i].dof_type = dof_int;
	copyHT(HT_parent_jointdisp0, frames[i].HT_parent_jd0);

	frames[i].iswheel = false;
	frames[i].issprocket = false;
	frames[i].isactuated = isactuated;
	frames[i].isfixed = isactuated;

	update_actframeinds();

}

void WmrModel::addWheelFrame(const std::string name, const std::string parent_name, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0, 
	const Real radius) {

	addJointFrame(name, parent_name, "RY", isactuated, HT_parent_jointdisp0);
	int i = nf-1;
	frames[i].iswheel = true;
	frames[i].rad = radius;

	update_wheelframeinds();
}

void WmrModel::addSprocketFrame(const std::string name, const std::string parent_name, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0, 
	const Real radius, const Real radius2, const Real len) {

	addJointFrame(name, parent_name, "RY", isactuated, HT_parent_jointdisp0);
	int i = nf-1;
	frames[i].issprocket = true;
	frames[i].rad = radius;
	frames[i].rad2 = radius2;
	frames[i].L = len;

	update_sprocketframeinds();
}

void WmrModel::setFrameMass(const int fi, const Real mass, const Vec3 center_of_mass, const Mat3 moment_of_inertia) {
	frames[fi].setMass(mass, center_of_mass, moment_of_inertia);
}

int WmrModel::nameToInd(const std::string name) const {
	int out = -1;
	for (int fi=0; fi<nf; fi++) {
		if (name.compare(frames[fi].name)==0) {
			out = fi;
			break;
		} 
	}
	return out;
};	


void inertiaBox(const Real m, const Real x, const Real y, const Real z, Mat3 I) {
	Vec3 Id;
	Id[0] = m*(y*y + z*z)/12;
	Id[1] = m*(x*x + z*z)/12;
	Id[2] = m*(x*x + y*y)/12;
	setDiagMat3(Id[0],Id[1],Id[2],I);
}

void inertiaCylinder(const Real m, const Real r, const Real h, const int axis, Mat3 I) {
	Vec3 Id;
	setcVec3(m*(3*r*r + h*h)/12,Id);
	Id[axis] = m*r*r/2;
	setDiagMat3(Id[0],Id[1],Id[2],I);
}

void KpKdToErpCfm( const Real kp, const Real kd, const Real dt, Real& erp, Real& cfm ) {
	erp = dt*kp/(dt*kp + kd);
	cfm = 1/(dt*kp + kd);
}
