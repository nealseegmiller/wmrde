#include <wmrde/surface/Surface.h>

/*
int surfacesHeight(const SurfaceVector& surfaces, const Vec3 pt, Real& height) {
	//OUTPUT
	//height:	height of surface si at point
	//(return)	si, index of surface for which height is highest

	Real height_; //temporary
	int loc;
	int surfidx;

	height = -REALMAX;
	for (int i=0; i<int(surfaces.size()); i++) { //loop over surfaces
		
		//assumes in bounds for at least one surface
		loc = surfaces[i]->surfaceHeight(pt,-1,height_);

		if ((loc >= 0) && (height_ > height)) {
			height = height_;
			surfidx=i;
		}
	}
	return surfidx;
}


int surfacesDz(const SurfaceVector& surfaces, const Vec3 pt, Real& dz, Vec3 normal) {
	//OUTPUT
	//dz:		distance of point from surface si
	//normal:	surface normal for surface si at point, can be null
	//(return)	si, index of surface for which dz is lowest

	Real dz_;
	Vec3 normal_;
	Real* ptr = 0; //can't set normal_ to zero
	if (normal != 0)
		ptr = normal_;
	int loc;
	int surfidx;

	
	dz = REALMAX;

	for (int i=0; i<int(surfaces.size()); i++) { //loop over surfaces

		//assumes in bounds for at least one surface
		loc = surfaces[i]->surfaceDz(pt,-1,dz_,ptr);

		//NaN less than comparison behavior is different in Debug and Release configurations?!
		//in Release configuration NaN < !NaN is true?!
		if ((loc >= 0) && (dz_ < dz)) {
			dz = dz_;
			if ( normal != 0 ) //if !null
				copyVec3(ptr,normal);
			surfidx = i;
		}

	}
	return surfidx;

}


*/





