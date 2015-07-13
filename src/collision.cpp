#include <wmrde/collision.h>

//update WheelContactGeom by discretizing wheel surface into points, selecting the point for which dz is minimized
void updateWheelContactGeomDiscretize(const SurfaceVector& surfaces, const HomogeneousTransform HT_wheel_to_world, const Real radius,  //input
	WheelContactGeom& contact) { //output

	const int MAXNP = 180+1; //max number of points in discretization
	
	Real angles[MAXNP]; //possible contact angles
	Vec3 pts_wheel[MAXNP]; //possible contact points along wheel surface (in wheel coords)
	Vec3 pts_world[MAXNP];
	Real deltazs[MAXNP]; //contact height errors for all possible pts
	int surfinds[MAXNP]; //surface indices for all possible pts
	int min_idx; //index of min value in deltazs
	int loc = -1; //TODO

	Real xproj = HT_wheel_to_world[2+COL0]; //dot(world z axis, wheel x axis)
	Real zproj = HT_wheel_to_world[2+COL2]; //dot(world z axis, wheel z axis)
	Real mid = atan2(xproj,zproj); // middle of range of possible contact angles

	Real rng = M_PI; //range
	int np=18+1; //number of points

	//discretize an arc of the wheel circumference to possible contact points
	linspace(mid-rng/2, mid+rng/2, np, angles);
	for (int i=0; i < np; i++) {
		contactAngleToPoint(radius, angles[i], pts_wheel[i]);
		applyHT(HT_wheel_to_world, pts_wheel[i], pts_world[i]);
	}

	//first step
	for (int i=0; i<np; i++) {
		surfinds[i] = surfacesDz(surfaces,pts_world[i],deltazs[i],0);
	}
	min_idx = findMin(np,deltazs);
	
	int si = surfinds[min_idx];

	if (1) {
		//second step

		mid = angles[min_idx];
		rng = M_PI/Real(np-1);
		np = 10+1;

		linspace(mid-rng/2, mid+rng/2, np, angles);
		for (int i=0; i<np; i++) {
			contactAngleToPoint(radius, angles[i], pts_wheel[i]);
			applyHT(HT_wheel_to_world,pts_wheel[i],pts_world[i]);
			surfaces[si]->surfaceDz(pts_world[i],loc,deltazs[i],0);
		}
		min_idx = findMin(np,deltazs);
	}

	contact.dz = deltazs[min_idx];
	contact.angle = angles[min_idx];

	Vec3 N_world, N_wheel; //normal vector
	surfaces[si]->surfaceNormal(pts_world[min_idx], loc, N_world);
	multMatTVec3(HT_wheel_to_world, N_world, N_wheel); //rotate to wheel coords

	HTContactToWheel(pts_wheel[min_idx], N_wheel, contact.HT_wheel ); //HT_contact_to_wheel
	composeHT(HT_wheel_to_world, contact.HT_wheel, contact.HT_world ); //HT_contact_to_world
	
}


//update WheelContactGeom by root finding. The dot product of a tangent vector to the wheel surface and the terrain normal is zero at the contact point.
//generally faster than the discretization method.
void updateWheelContactGeomRoot(const SurfaceVector& surfaces, const HomogeneousTransform HT_wheel_to_world, const Real radius, 
	WheelContactGeom& contact) { //output

	//http://en.wikipedia.org/wiki/Brent's_method
	
	Real xproj = HT_wheel_to_world[2+COL0]; //dot(world z axis, wheel x axis)
	Real zproj = HT_wheel_to_world[2+COL2]; //dot(world z axis, wheel z axis)
	Real mid = atan2(xproj,zproj); // middle of range of possible contact angles

	Real rng = (2.0/3.0)*M_PI; //range

	Real a = mid - rng/2;
	Real b = mid + rng/2;

	Real tolx = DEGTORAD(.5);
	Real tolfx = 1e-2;
	
	
	Real dz;
	Vec3 pt_wheel; //contact point (in wheel coords)
	Vec3 N_world; //terrain normal (in world coords)

	//lambda closure, requires C++11
	//alternative is functor
	//using auto is faster than std::function
	auto ddz_dca = [&] (const Real angle) -> Real {
		
		//captures surfaces, HT_wheel_to_world, radius
		//also captures dz, pt_wheel, normal ?

		Vec3 pt_world;
		Vec3 tangent; //wheel tangent (in world coords)

		contactAngleToPoint(radius, angle, pt_wheel);
		applyHT(HT_wheel_to_world, pt_wheel, pt_world);

		surfacesDz(surfaces,pt_world,dz,N_world);

		Real s = sin(angle);
		Real c = cos(angle);
		for (int i=0; i<3; i++) 
			tangent[i] = - c*HT_wheel_to_world[i] + s*HT_wheel_to_world[i+COL2];

		Real f = dotVec3(tangent,N_world);
		return f;
	};

	Real contact_angle;
	//uncomment *one* of the following:
	//contact_angle = findRootBisection(a,b,ddz_dca,tolx,tolfx);
	contact_angle = findRootBrents(a,b,ddz_dca,tolx,tolfx); 

	contact.dz = dz;
	contact.angle = contact_angle;
	
	Vec3 N_wheel;
	multMatTVec3(HT_wheel_to_world, N_world, N_wheel); //rotate to wheel coords

	HTContactToWheel(pt_wheel, N_wheel, contact.HT_wheel ); //HT_contact_to_wheel
	composeHT(HT_wheel_to_world, contact.HT_wheel, contact.HT_world ); //HT_contact_to_world
}


//compute transform from contact to wheel coords given contact point & terrain normal (in wheel coords)
void HTContactToWheel(const Vec3 pt, const Vec3 N, HomogeneousTransform HT_contact_to_wheel) {

	copyVec3(N,HT_contact_to_wheel+COL2); //z axis = terrain normal

	Vec3 wheel_yaxis;
	setVec3(0,1,0,wheel_yaxis);
	crossVec3(wheel_yaxis,HT_contact_to_wheel+COL2,HT_contact_to_wheel); //x axis = cross(wheel frame y axis, z axis)
	normalizeVec3(HT_contact_to_wheel);

	crossVec3(HT_contact_to_wheel+COL2,HT_contact_to_wheel,HT_contact_to_wheel+COL1); //y axis =  cross(z axis, x axis)
	copyVec3(pt,HT_contact_to_wheel+COL3); //origin is contact point
}

void whichPointsInContactWheel(const int min_npic, const int nw, WheelContactGeom* contacts) {
	if (min_npic == 0) {
		//dynamic simulation
		for (int wno=0; wno < nw; wno++)
			contacts[wno].incontact = contacts[wno].dz < 0;
	} else {
		//kinematic simulation
		//TODO
		for (int wno=0; wno < nw; wno++) //wheel number
			contacts[wno].incontact = true;
	}
}

void whichPointsInContactTrack(const int min_npic, const int nt, TrackContactGeom* contacts) {
	if (min_npic == 0) {
		//dynamic simulation
		for (int tno=0; tno < nt; tno++) {
			int np = contacts[tno].get_np(); //number of points
			for (int pno=0; pno < np; pno++)
				contacts[tno].incontact[pno] = contacts[tno].dz[pno] < 0;
		}
	} else {
		//kinematic simulation
		//TODO, not in contact if dz is NaN!
		for (int tno=0; tno < nt; tno++) //track number
			setVec(contacts[tno].get_np(), true, contacts[tno].incontact);
	}	
}

void whichPointsInContact(const int min_npic, const int np, const Vec3 pts_world[], const Real dz[], int incontact[]) {
	//TODO
	//determine stability based on eigenvalues of covariance of points
}

//for tracks
void initTrackContactGeom(const Real rad, const Real rad2, const Real L, TrackContactGeom& contact) {
	const int npflat = 3;
	const bool sides[4] = {true, false, false, false};
	//const bool sides[4] = {true, true, true, true};

	int np;
	np = HTContactToTrack(rad,rad2,L,npflat,sides,contact.HT_track);
	contact.set_np(np);	
}

int HTContactToTrack(const Real rad1, const Real rad2, const Real L, const int npflat, const bool sides[], HomogeneousTransform HT_track[]) {
	int np = 0; //number of contact points

	Real ca = -asin((rad1-rad2)/L); //contact angle
	Real Lflat = L*cos(ca); //length of flat segment
	Real dx = Lflat/(Real(npflat)-1);

	if (sides[0]) { //lower flat
		
		Vec3 cp; //contact point
		contactAngleToPoint(rad1, ca, cp);

		Vec3 N = {0,0,0}; //normal
		addmVec3(N,-1.0,cp,N);
		normalizeVec3(N);

		HomogeneousTransform HT; //temporary
		HTContactToWheel(cp, N, HT);
		Real* xaxis = HT; //1st col

		for (int i=0; i<npflat; i++) {
			if (i>0)
				addmVec3(cp,dx,xaxis,cp);

			if (i == 0 && sides[3]) //avoid duplicate points
				continue;

			assert(np < ContactGeom::MAXNP);
			HTContactToWheel(cp, N, HT_track[np]);
			np++;
		}
	}

	Real rca = REALSIGN(L)*M_PI + 2*ca; //range of contact angles for sprocket

	if (sides[1]) { //front sprocket
		int npsprocket = int(round(rca*rad2/dx)) + 1;
		Real dca = rca/(Real(npsprocket) - 1);
		Real angle = ca;

		for (int i=0; i<npsprocket; i++) {
			if (i>0)
				angle -= dca;

			if (i == 0 && sides[0]) //avoid duplicate points
				continue;

			Vec3 cp; //contact point
			contactAngleToPoint(rad2,angle,cp);
			cp[0] += L;

			Vec3 N = {L,0,0}; //normal
			addmVec3(N,-1.0,cp,N);
			normalizeVec3(N);

			assert(np < ContactGeom::MAXNP);
			HTContactToWheel(cp,N,HT_track[np]);
			np++;
		}
	}

	if (sides[2]) { //upper flat
		
		Vec3 cp; //contact point
		contactAngleToPoint(rad2, ca-rca, cp);
		cp[0] += L;

		Vec3 N = {L,0,0}; //normal
		addmVec3(N,-1.0,cp,N);
		normalizeVec3(N);

		HomogeneousTransform HT; //temporary
		HTContactToWheel(cp, N, HT);
		Real* xaxis = HT; //1st col

		for (int i=0; i<npflat; i++) {
			if (i>0)
				addmVec3(cp,dx,xaxis,cp);

			if (i == 0 && sides[1]) //avoid duplicate points
				continue;

			assert(np < ContactGeom::MAXNP);
			HTContactToWheel(cp, N, HT_track[np]);
			np++;
		}
	}

	rca = REALSIGN(L)*M_PI - 2*ca;

	if (sides[3]) { //rear sprocket
		int npsprocket = int(round(rca*rad1/dx)) + 1;
		Real dca = rca/(Real(npsprocket) - 1);
		Real angle = ca+rca;

		for (int i=0; i<npsprocket; i++) {
			if (i>0)
				angle -= dca;

			if (i == 0 && sides[2]) //avoid duplicate points
				continue;

			Vec3 cp; //contact point
			contactAngleToPoint(rad1,angle,cp);

			Vec3 N = {0,0,0}; //normal
			addmVec3(N,-1.0,cp,N);
			normalizeVec3(N);

			assert(np < ContactGeom::MAXNP);
			HTContactToWheel(cp,N,HT_track[np]);
			np++;
		}
	}

	return np;
}


void updateTrackContactGeom(const SurfaceVector& surfaces, const HomogeneousTransform HT_track_to_world, //input
	TrackContactGeom& contact) { //output

	int np = contact.get_np();

	for (int pno=0; pno<np; pno++) {
		//HT_contact_to_world
		composeHT(HT_track_to_world, contact.HT_track[pno], contact.HT_world[pno]);

		Real* cp_world = contact.HT_world[pno] + COL3; //contact point in world coords
		Real Nz_world = contact.HT_world[pno][COL2+2]; //z component of normal vector in world coords

		Real zsurf;
		surfacesHeight(surfaces, cp_world, zsurf);

		Real dh = cp_world[2] - zsurf; //delta height
		if (Nz_world > 0) {
			contact.dz[pno] = dh * Nz_world;
		} else {
			contact.dz[pno] = REALNAN; //don't use these
		}
	}
}
