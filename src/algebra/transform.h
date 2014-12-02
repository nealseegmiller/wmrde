//transform.h
//functions for coordinate transformation in R3
//Homogeneous transform is 4x4 matrix: 3x4 matrix = [R|t] on top, bottom row = [0 0 0 1]
//only need to store 12 unique elements of [R|t]
//extra elements for 16 byte alignment
//can convert between macros and inline functions based on profiling

#ifndef _WMRSIM_TRANSFORM_H_
#define _WMRSIM_TRANSFORM_H_

#include <algebra/linalg3.h>

typedef Real HomogeneousTransform[4*SIZEVEC3];

//HT indices:
//(column-major order)
//0,4,8,  12
//1,5,9,  13
//2,6,10, 14
//3,7,11, 15 //unused, for 16 byte alignment only

//r11,r12,r13, t1
//r21,r22,r23, t2
//r31,r32,r33, t3

#define COL3 12

//if HT is pointer to first element of Homogeneous Transform
//HT is also pointer to Rotation matrix
//HT+COL3 is pointer to translation vector

//orientation vectors
#define SIZEEULER 3
#define SIZEQUAT 4

typedef Real VecEuler[SIZEEULER]; //[roll pitch yaw]
typedef Real MatEuler[3*SIZEEULER]; //for 3x3 transforms between Euler rates and angular velocity
typedef Real VecQuat[SIZEQUAT]; //[w x y z]
typedef Real MatQuat[3*SIZEQUAT]; //for 3x4 or 4x3 transforms between quaternion rates and angular velocity

#define setIdentityHT(Dest) do { \
	setIdentityMat3(Dest); \
	setcVec3(0.0,Dest + COL3); \
} while(0)

//inline void copyHT(const HomogeneousTransform Source, HomogeneousTransform Dest) {
#define copyHT(Source,Dest) do { \
	copyMat3(Source,Dest); \
	copyVec3(Source+COL3,Dest+COL3); \
} while(0)


//ROTATION FUNCTIONS
//inline void Rotx(const Real a, Mat3 R) {
#define Rotx(a,R) do { \
	R[0] = 1; R[0+COL1] = 0;			R[0+COL2] = 0; \
	R[1] = 0; R[1+COL1] = cos(a);		R[1+COL2] =-sin(a); \
	R[2] = 0; R[2+COL1] =-R[1+COL2];	R[2+COL2] = R[1+COL1]; \
} while(0)

//inline void Roty(const Real a, Mat3 R) {
#define Roty(a,R) do { \
	R[0] = cos(a);		R[0+COL1] = 0;	R[0+COL2] = sin(a); \
	R[1] = 0;			R[1+COL1] = 1;	R[1+COL2] = 0; \
	R[2] =-R[0+COL2];	R[2+COL1] = 0;	R[2+COL2] = R[0]; \
} while(0)


//inline void Rotz(const Real a, Mat3 R) {
#define Rotz(a,R) do { \
	R[0] = cos(a);		R[0+COL1] =-sin(a);	R[0+COL2] = 0; \
	R[1] =-R[0+COL1];	R[1+COL1] = R[0];	R[1+COL2] = 0; \
	R[2] = 0;			R[2+COL1] = 0;		R[2+COL2] = 1; \
} while(0)


//Rotation about axis and angle
//i:	axis index
//a:	angle
//R:	Rotation matrix
//inline void RotInd(const int i, const Real angle, Mat3 R) {
#define RotInd(i,angle,R) do { \
	if (i==0) \
		Rotx(angle,R); \
	else if (i==1) \
		Roty(angle,R); \
	else if (i==2) \
		Rotz(angle,R); \
} while(0)


//to invert HT:
//R_inv = R'
//t_inv = -R'*t
//inline void invertHT(const HomogeneousTransform T, HomogeneousTransform invT) {
#define invertHT(T,invT) do { \
	copyTMat3(T,invT); \
	multMatTVec3(T,T+COL3,invT+COL3); \
	mulcVec3(invT+COL3,-1,invT+COL3); \
} while(0)


//COMPOSE FUNCTIONS

//to compose HT C=A*B:
//Rc = Ra*Rb
//tc = Ra*tb + ta
inline void composeHT(const HomogeneousTransform A, const HomogeneousTransform B, HomogeneousTransform C) {
//#define composeHT(A,B,C) do { 
	multMatMat3(A,B,C); \
	multMatVec3(A,B+COL3,C+COL3); \
	addVec3(A+COL3,C+COL3,C+COL3); \
} //while(0)


//to compose HT C=inv(A)*B:
//Rc = Ra'*Rb
//tc = Ra'*tb - Ra'*ta
//inline void composeInvHT(const HomogeneousTransform A, const HomogeneousTransform B, HomogeneousTransform C) {
#define composeInvHT(A,B,C) do { \
	multMatTVec3(A,B+COL3,C+COL3); \
	multMatTVec3(A,A+COL3,C); /*temporary*/\
	addmVec3(C+COL3,-1,C,C+COL3); \
	multMatTMat3(A,B,C); \
} while(0)

//APPLY FUNCTIONS

//transform coordinates of a point
//q = R*p + t
//inline void applyHT(const HomogeneousTransform T, const Vec3 p, Vec3 q) {
#define applyHT(T,p,q) do { \
	multMatVec3(T,p,q); \
	addVec3(q,T+COL3,q); \
} while(0)

//use inverse of homogeneous transform
//q = R'*p - R'*t = R'*(p-t)
inline void applyInvHT(const HomogeneousTransform T, const Vec3 p, Vec3 q) {
	Vec3 tmp; //temporary
	addmVec3(p,-1,T+COL3,tmp);
	multMatTVec3(T,tmp,q);
}

//PRINT FUNCTION
inline void printHT(const HomogeneousTransform T, int precision, int width) {
	Real T_[3*4];
	copyVec3ToArray(4,(const Vec3*) T,T_);
	printMatReal(3,4,T_,precision,width);
}


//ORIENTATION FUNCTIONS
//functions for both Euler angles and quaternions should be available regardless of which is used for simulation

inline void setEuler(const Real rol, const Real pit, const Real yaw, VecEuler e) {
	e[0] = rol;
	e[1] = pit;
	e[2] = yaw;
}
inline void copyEuler(const VecEuler Source, VecEuler Dest) {
	for (int i=0; i<SIZEEULER; i++)
		Dest[i] = Source[i];
}

inline void setQuat( const Real w, const Real x, const Real y, const Real z, VecQuat q) {
	q[0]=w;
	q[1]=x;
	q[2]=y;
	q[3]=z;
}
inline void copyQuat( const VecQuat Source, VecQuat Dest) {
	for (int i=0; i<SIZEQUAT; i++)
		Dest[i] = Source[i];
}

inline void normalizeQuat( VecQuat q) {
	Real n = normVec3(q+1);
	for (int i=0; i<SIZEQUAT; i++)
		q[i] /= n;
}


void eulerToRot(const VecEuler euler, Mat3 R); //called once per time step
void rotToEuler(const Mat3 R, VecEuler euler); //called infrequently
void velToEulerrate(const VecEuler euler, const Vec3 vel, VecEuler eulerrate, MatEuler T); //called once per time step
void eulerrateToVel(const VecEuler euler, const VecEuler eulerrate, Vec3 vel, MatEuler T); //called infrequently

//TODO
void quatToRot(const VecQuat quat, Mat3 R);
void rotToQuat(const Mat3 R, VecQuat quat);
void velToQuatrate(const VecQuat quat, const Vec3 vel, VecQuat quatrate, MatQuat T);
void quatrateToVel(const VecQuat quat, const VecQuat quatrate, Vec3 vel, MatQuat T);

void eulerToQuat(const VecEuler euler, VecQuat quat);
void quatToEuler(const VecQuat quat, VecEuler euler);

//PRINT FUNCTIONS
inline void printEuler(const VecEuler euler, int precision, int width) { 
	printMatReal(SIZEEULER,1,euler,precision,width); 
}
inline void printQuat(const VecQuat quat, int precision, int width) { 
	printMatReal(SIZEQUAT,1,quat,precision,width); 
}


//generic orient types to easily switch between Euler angles and quaternions
#if WMRSIM_USE_QUATERNION
typedef VecQuat VecOrient;
typedef MatQuat MatOrient;
#define copyOrient(source,dest) copyQuat(source,dest)
#define orientToRot(orient,R) quatToRot(orient,R)
#define rotToOrient(R,orient) rotToQuat(R,orient)
#define velToOrientrate(orient,vel,orientrate,T) velToQuatrate(orient,vel,orientrate,T)
#define orientrateToVel(orient,orientrate,vel,T) quatrateToVel(orient,orientrate,vel,T)
#define printOrient(orient,precision,width) printQuat(orient,precision,width)
#else
typedef VecEuler VecOrient;
typedef MatEuler MatOrient;
#define copyOrient(source,dest) copyEuler(source,dest)
#define orientToRot(orient,R) eulerToRot(orient,R)
#define rotToOrient(R,orient) rotToEuler(R,orient)
#define velToOrientrate(orient,vel,orientrate,T) velToEulerrate(orient,vel,orientrate,T)
#define orientrateToVel(orient,orientrate,vel,T) eulerrateToVel(orient,orientrate,vel,T)
#define printOrient(orient,precision,width) printEuler(orient,precision,width)
#endif

//convert pose (orientation and translation vectors) to homogeneous transform
//inline void poseToHT(const VecOrient orient, const Vec3 t, HomogeneousTransform T) {
#define poseToHT(orient,t,T) do { \
	orientToRot(orient,T); \
	copyVec3(t,T+COL3); \
} while(0)

//convert homogeneous transform to pose
//inline void HTToPose(const HomogeneousTransform T, VecOrient orient, Vec3 t) {
#define HTToPose(T,orient,t) do { \
	rotToOrient(T,orient); \
	copyVec3(T+COL3,t); \
} while(0)



#endif //_WMRSIM_TRANSFORM_H_