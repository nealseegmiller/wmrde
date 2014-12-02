#include <surface/TriMeshSurf.h>

TriMeshSurf::TriMeshSurf() {
	nv = 0;
	nt = 0;

	xmin = REALMAX;
	xmax = -REALMAX;
	ymin = REALMAX;
	ymax = -REALMAX;
}

void TriMeshSurf::add_vertex(const Vec3 V) {
	copyVec3(V,vertices[nv]);
	nv++;

	//bounding box
	if (V[0] < xmin)
		xmin = V[0];
	if (V[0] > xmax)
		xmax = V[0];
	if (V[1] < ymin)
		ymin = V[1];
	if (V[1] > ymax)
		ymax = V[1];
}

void TriMeshSurf::add_triangle(const int i0, const int i1, const int i2) {
	indices[S2I(0,nt,3)] = i0;
	indices[S2I(1,nt,3)] = i1;
	indices[S2I(2,nt,3)] = i2;

	//set pec
	Vec3 v01,v02,N; //intermediate
	addmVec3(vertices[i1],-1,vertices[i0],v01); //vector from 0 to 1
	addmVec3(vertices[i2],-1,vertices[i0],v02); //vector from 0 to 2
	crossVec3(v01,v02,N); //normal
	normalizeVec3(N);
	
	Real* pec_ = pec + S2I(0,nt,4); //pointer to col
	copyVec3(N,pec_);
	pec_[3] = -dotVec3(pec_,vertices[i0]);

	//set ineq
	get_ineq(vertices[i0],vertices[i1],M01[nt]);
	get_ineq(vertices[i1],vertices[i2],M12[nt]);
	get_ineq(vertices[i2],vertices[i0],M20[nt]);

	nt++;
}

//add a quad, comprising two triangles
void TriMeshSurf::add_quad(const int i0, const int i1, const int i2, const int i3) {
	add_triangle(i0,i1,i2);
	add_triangle(i0,i2,i3);
}

//get location
//(return)		loc, index of triangle, -1 if out of bounds
int TriMeshSurf::get_loc(const Vec3 pt) {
	
	int loc = -1;
	if (pt[0] >= xmin && pt[0] <= xmax && 
		pt[1] >= ymin && pt[1] <= ymax) {

		Vec3 pt_; //temporary
		pt_[0] = pt[0]; 
		pt_[1] = pt[1]; 
		pt_[2] = 1.0;

		//test if in bounds for each triangle
		for (int i=0; i<nt; i++) {
			if (dotVec3(M01[i],pt_) <= 0 && 
				dotVec3(M12[i],pt_) <= 0 && 
				dotVec3(M20[i],pt_) <= 0) {
				loc = i;
				break;
			}
		}
	}
	return loc;
}

//get inequality
//lec:		line equation coefficients for the line through the vertices, [a b c]*[x y 1]^T = 0
void TriMeshSurf::get_ineq(const Vec3 v0, const Vec3 v1, Vec3 lec) {

	Real a,b,c;
	a = v1[1] - v0[1]; //dy
	b =-(v1[0] - v0[0]); //-dx

	//normalize [a b]
	Real n = sqrt(a*a + b*b);
	a = a/n; b = b/n;

	c = -(a*v0[0] + b*v0[1]);

	lec[0] = a;
	lec[1] = b;
	lec[2] = c;
}



int TriMeshSurf::surfaceHeight(const Vec3 pt, int loc, Real& height) {
	if (loc < 0)
		loc = get_loc(pt);

	if (loc >= 0) {
		Real* pec_ = pec + S2I(0,loc,4); //pointer to col

		//z = -1/c*[a b d]*[x y 1]'
		height = -(pec_[0]*pt[0] + pec_[1]*pt[1] + pec_[3])/pec_[2];
	} else {
		height = REALNAN;
	}

	return loc;
}

int TriMeshSurf::surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal) {
	if (loc < 0)
		loc = get_loc(pt);

	if (loc >= 0) {
		Real* pec_ = pec + S2I(0,loc,4); //pointer to col

		//dz = [a b c d]*[x y z 1]'
		dz = pec_[0]*pt[0] + pec_[1]*pt[1] + pec_[2]*pt[2] + pec_[3];

		if (normal != 0)
			copyVec3(pec_,normal);
	} else {
		dz = REALNAN;
	}

	return loc;
}

int TriMeshSurf::surfaceNormal(const Vec3 pt, int loc, Vec3 normal) {
	if (loc < 0)
		loc = get_loc(pt);

	if (loc >= 0)
		copyVec3(pec+S2I(0,loc,4),normal);

	return loc;
}