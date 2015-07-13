#include <wmrde/surface/PlaneSurf.h>

//overloaded constructor
PlaneSurf::PlaneSurf(const Real PlaneEqnCoeff[4]) {
	Real norm = normVec3(PlaneEqnCoeff); //first three elements of pec should be unit normal vector
	for (int i=0; i<4; i++)
		pec[i] = PlaneEqnCoeff[i]/norm;
}

PlaneSurf::PlaneSurf(const Vec3 normal, const Vec3 point) {
	copyVec3(normal,pec);
	normalizeVec3(pec);
	pec[3] = -dotVec3(pec,point);
}



//define the pure virtual equations required by the Abstract super class
int PlaneSurf::surfaceHeight(const Vec3 pt, int loc, Real& height) {
	//z = -1/c*[a b d]*[x y 1]'
	height = -(pec[0]*pt[0] + pec[1]*pt[1] + pec[3])/pec[2];

	return 0;
}

int PlaneSurf::surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal) {

	//dz = [a b c d]*[x y z 1]'
	dz = pec[0]*pt[0] + pec[1]*pt[1] + pec[2]*pt[2] + pec[3];

	if (normal != 0)
		copyVec3(pec,normal);

	return 0;
}

int PlaneSurf::surfaceNormal(const Vec3 pt, int loc, Vec3 normal) {
	copyVec3(pec,normal);
	return 0;
}




