//TriMeshSurf.h
//The TriMeshSurf class is derived from the Surface class
//This class specifies a triangular mesh

#ifndef _WMRSIM_TRIMESHSURF_H_
#define _WMRSIM_TRIMESHSURF_H_

#include <surface/Surface.h>
#include <algebra/transform.h>

class TriMeshSurf : public Surface {
public:
	static const int MAXNT = 20; //max number of triangles

private:
	int nv;
	int nt;

	Vec3 vertices[3*MAXNT];
	int indices[3*MAXNT];	//3 x nt, 3 indices per triangle
	Real pec[4*MAXNT];		//4 x nt, plane equation coefficients

	//inequality matrices, one for each edge of triangle: 0-1, 1-2, 2-0
	//used to test if point is in bounds of triangle
	Vec3 M01[MAXNT];
	Vec3 M12[MAXNT];
	Vec3 M20[MAXNT];

	Real xmin,xmax,ymin,ymax;	//bounding box

	int get_loc(const Vec3 pt);
	void get_ineq(const Vec3 v0, const Vec3 v1, Vec3 lec);

public:
	TriMeshSurf();

	int get_nv() const { return nv; }
	int get_nt() const { return nt; }
	const Vec3* get_vertices() const { return vertices; }
	const int* get_indices() const { return indices; }
	const Real* get_pec() const { return pec; }

	void add_vertex(const Vec3 V);
	//specify indices in counter-clockwise order!
	void add_triangle(const int i0, const int i1, const int i2); 
	void add_quad(const int i0, const int i1, const int i2, const int i3);

	int surfaceHeight(const Vec3 pt, int loc, Real& height);
	int surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal);
	int surfaceNormal(const Vec3 pt, int loc, Vec3 normal);

};


#endif