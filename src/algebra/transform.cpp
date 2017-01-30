#include <wmrde/algebra/transform.h>
#include <wmrde/util/common_util.h>

namespace wmrde
{

//euler = [rol pit yaw]
//R = Rotz(yaw)*Roty(pit)*Rotx(rol)
void eulerToRot(const VecEuler euler, Mat3 R) {

	Real sR=sin(euler[0]); Real cR=cos(euler[0]); //rol
	Real sP=sin(euler[1]); Real cP=cos(euler[1]); //pit
	Real sY=sin(euler[2]); Real cY=cos(euler[2]); //yaw
	R[0]=cP*cY;	R[0+COL1]=cY*sP*sR - cR*sY;	R[0+COL2]=sR*sY + cR*cY*sP;
	R[1]=cP*sY;	R[1+COL1]=cR*cY + sP*sR*sY;	R[1+COL2]=cR*sP*sY - cY*sR;
	R[2]=-sP;	R[2+COL1]=cP*sR;			R[2+COL2]=cP*cR;
}

void rotToEuler(const Mat3 R, VecEuler euler) {
	euler[2] = atan2(R[1], R[0]); //yaw
	euler[0] = atan2(R[2+COL1], R[2+COL2]); //rol
	Real sR=sin(euler[0]);
	Real cR=cos(euler[0]);
	Real cP;
	if (fabs(cR) > fabs(sR)) {
		cP = R[2+COL2]/cR; //singularity at rol = +/- pi/2
	} else {
		cP = R[2+COL1]/sR; //singularity at rol = 0,+/-pi
	}
	euler[1] = atan2(-R[2], cP); //pit
}

//eulerrate=T*vel
void velToEulerrate(const VecEuler euler, const Vec3 vel, VecEuler eulerrate, MatEuler T) {
	//singularity at pitch = +/- pi/2, rol & yaw axes align
	//problems if single precision

	/*
	Real sR = sin(euler[0]);
	Real cR = cos(euler[0]);
	Real cP = cos(euler[1]);
	Real tP = tan(euler[1]);
	*/
	
	double sR = sin(double(euler[0]));
	double cR = cos(double(euler[0]));
	double cP = cos(double(euler[1]));
	double tP = tan(double(euler[1]));

	//may be null parameters
	if (T) {
		const int col1 = SIZEEULER;
		const int col2 = 2*SIZEEULER;

		T[0]=1;	T[0+col1]=sR*tP; T[0+col2]=cR*tP;
		T[1]=0;	T[1+col1]=cR;	 T[1+col2]=-sR;
		T[2]=0; T[2+col1]=sR/cP; T[2+col2]=cR/cP;
	}
	if (vel && eulerrate) {
		eulerrate[0] = vel[0] + sR*tP*vel[1] + cR*tP*vel[2];
		eulerrate[1] = cR*vel[1] - sR*vel[2];
		eulerrate[2] = sR/cP*vel[1] + cR/cP*vel[2];
	}
}

//vel=T*eulerrate
void eulerrateToVel(const VecEuler euler, const VecEuler eulerrate, Vec3 vel, MatEuler T) {
	Real sR = sin(euler[0]);
	Real cR = cos(euler[0]);
	Real sP = sin(euler[1]);
	Real cP = cos(euler[1]);
	//may be null parameters
	if (T) {
		const int col1 = SIZEEULER;
		const int col2 = 2*SIZEEULER;

		T[0]=1;	T[0+col1]=0;	T[0+col2]=-sP;
		T[1]=0;	T[1+col1]=cR;	T[1+col2]=sR*cP;
		T[2]=0; T[2+col1]=-sR;	T[2+col2]=cR*cP;
	}
	if (eulerrate && vel) {
		vel[0] = eulerrate[0] - sP*eulerrate[2];
		vel[1] = cR*eulerrate[1] + sR*cP*eulerrate[2];
		vel[2] = -sR*eulerrate[1] + cR*cP*eulerrate[2];
	}
}

//TODO
void quatToRot(const VecQuat quat, Mat3 R) {
	//TODO, matrix multiplication method?
	//Real m;
	//NORMALIZEQUAT(quat,m)

	Real w=quat[0];
	Real x=quat[1];
	Real y=quat[2];
	Real z=quat[3];

	Real wx = w*x; 
	Real wy = w*y; 
	Real wz = w*z;
	Real xx = x*x; 
	Real xy = x*y; 
	Real xz = x*z;
	Real yy = y*y; 
	Real yz = y*z; 
	Real zz = z*z;

	R[0] = 1.0-2.0*(yy+zz); R[0+COL1] = 2.0*(xy-wz);	 R[0+COL2] = 2.0*(xz+wy);
	R[1] = 2.0*(xy+wz);		R[1+COL1] = 1.0-2.0*(xx+zz); R[1+COL2] = 2.0*(yz-wx);
	R[2] = 2.0*(xz-wy);		R[2+COL1] = 2.0*(yz+wx);	 R[2+COL2] = 1.0-2.0*(xx+yy);

}
void rotToQuat(const Mat3 R, VecQuat quat) {
	Real Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz;
	Qxx = R[0]; Qxy = R[0+COL1]; Qxz = R[0+COL2];
	Qyx = R[1]; Qyy = R[1+COL1]; Qyz = R[1+COL2];
	Qzx = R[2]; Qzy = R[2+COL1]; Qzz = R[2+COL2];

	Real t = Qxx+Qyy+Qzz;
	Real r = sqrt(1+t);
	Real w = 0.5*r;

	Real x,y,z;
	if (t > 0) {
		Real s = 0.5/r;
		x = (Qzy-Qyz)*s;
		y = (Qxz-Qzx)*s;
		z = (Qyx-Qxy)*s;
	} else {
		//avoid singularity
		//use std::copysign if C++11
		x = REALSIGN(Qzy-Qyz)*0.5*sqrt(1+Qxx-Qyy-Qzz);
		y = REALSIGN(Qxz-Qzx)*0.5*sqrt(1-Qxx+Qyy-Qzz);
		z = REALSIGN(Qyx-Qxy)*0.5*sqrt(1-Qxx-Qyy+Qzz);
	}
	quat[0] = w; quat[1] = x; quat[2] = y; quat[3] = z;
}
//quatrate=T*vel
void velToQuatrate(const VecQuat quat, const Vec3 vel, VecQuat quatrate, MatQuat T) {
	//TODO
}
//vel=T*quatrate
void quatrateToVel(const VecQuat quat, const VecQuat quatrate, Vec3 vel, MatQuat T) {
	//TODO
}

void eulerToQuat(const VecEuler euler, VecQuat quat) {
	//TODO, a more direct way
	Mat3 R;
	eulerToRot(euler,R);
	rotToQuat(R,quat);
}
void quatToEuler(const VecQuat quat, VecEuler euler) {
	//TODO, a more direct way
	Mat3 R;
	quatToRot(quat,R);
	rotToEuler(R,euler);
}

} //namespace
