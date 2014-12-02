#include <surface/GridSurf.h>

int GridSurf::readfile(const std::string FileName ) {
	// declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
	std::ifstream file ( FileName.c_str() );

	std::string value;
	if ( file.is_open() ) {
		for (int i=0; i<2; i++) {
			std::getline (file, value, ','); lowerlim[i] = (Real) atof(value.c_str());
			std::getline (file, value, ','); upperlim[i] = (Real) atof(value.c_str());
			std::getline (file, value, ','); n[i] = atoi(value.c_str());

			assert(upperlim[i] > lowerlim[i]);
			d[i] = (upperlim[i]-lowerlim[i])/(n[i]-1);
		}

		int i=0;
		while ( file.good() ) {
			assert( i < buffer ); //check if buffer exceeded
			std::getline ( file, value, ',' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
			Z[i] = (Real) atof(value.c_str());
			i++;
		}
		file.close();
	
		assert(i == (n[0]*n[1])+1); //check if correct number of elements in Z
		
		calcGradients();
		return 1;
	} else {
		return 0;
	}

}

void GridSurf::calcGradients() {
	//called by constructor
	for (int i=0; i<n[0]; i++) {
		for (int j=0; j<n[1]; j++) {

			int im,ip,jm,jp;
			im=std::max(i-1,0);
			ip=std::min(i+1,n[0]-1);
			jm=std::max(j-1,0);
			jp=std::min(j+1,n[1]-1);

			DZDX[S2I2(i,j,n[0])] = (Z[S2I2(ip,j,n[0])] - Z[S2I2(im,j,n[0])]) / ((ip-im)*d[0]);
			DZDY[S2I2(i,j,n[0])] = (Z[S2I2(i,jp,n[0])] - Z[S2I2(i,jm,n[0])]) / ((jp-jm)*d[1]);
		}
	}
}

int GridSurf::get_loc(const Vec3 pt) {
	//check if in bounds
	if (pt[0] < lowerlim[0] || pt[0] > upperlim[0] || pt[1] < lowerlim[1] || pt[1] > upperlim[1]) 
		return -1; 
	else 
		return 0;
}

int GridSurf::surfaceHeight(const Vec3 pt, int loc, Real& height) {
	if (loc < 0)
		loc = get_loc(pt);

	if (loc >= 0) {
		
		int I[2];
		Real xn[2],f[4];
		bilinearInterp_xn(lowerlim, d, pt, I, xn);
		bilinearInterp_f(Z, I, n, f);
		height = bilinearInterp(f,xn);
		
	} else {
		height = REALNAN;
	}
	return loc;
}

int GridSurf::surfaceDz(const Vec3 pt, int loc, Real& dz, Vec3 normal) {
	if (loc < 0)
		loc = get_loc(pt);

	if (loc >= 0) {

		//dot product of delta height (along world z axis) and surface normal
		Real dzdx,dzdy;
		
		int I[2];
		Real xn[2],f[4];
		bilinearInterp_xn(lowerlim, d, pt, I, xn);

		bilinearInterp_f(Z, I, n, f);
		dz = bilinearInterp(f,xn);
		bilinearInterp_f(DZDX, I, n, f);
		dzdx = bilinearInterp(f,xn);
		bilinearInterp_f(DZDY, I, n, f);
		dzdy = bilinearInterp(f,xn);
		
		Vec3 N;
		N[0] = -dzdx;
		N[1] = -dzdy;
		N[2] = 1.0;
		normalizeVec3(N);

		dz = pt[2] - dz;
		dz = dz*N[2];

		if (normal != 0)
			copyVec3(N,normal);

	} else {
		dz = REALNAN;
	}
	return loc;
}

int GridSurf::surfaceNormal(const Vec3 pt, int loc, Vec3 normal) {
	if (loc < 0)
		loc = get_loc(pt);

	if (loc >= 0) {
		Real dzdx,dzdy;
		
		int I[2];
		Real xn[2],f[4];
		bilinearInterp_xn(lowerlim, d, pt, I, xn);

		bilinearInterp_f(DZDX, I, n, f);
		dzdx = bilinearInterp(f,xn);
		bilinearInterp_f(DZDY, I, n, f);
		dzdy = bilinearInterp(f,xn);
		

		normal[0] = -dzdx;
		normal[1] = -dzdy;
		normal[2] = 1.0;
		normalizeVec3(normal);
	}
	return loc;
}