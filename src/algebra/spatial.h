//spatial.h
//functions for spatial vector algebra
//Reference: Roy Featherstone, Rigid Body Dynamics Algorithms
//http://royfeatherstone.org/spatial/index.html#spatial-software
//can convert between macros and inline functions based on profiling

#ifndef _WMRSIM_SPATIAL_H_
#define _WMRSIM_SPATIAL_H_

#include <algebra/transform.h>
#include <algebra/matrix.h>

typedef Real Vec6b[2*SIZEVEC3];
typedef Real Mat6b[4*SIZEMAT3];

//spatial vectors are 6x1 (two 3x1 blocks)
//0
//1
//2
//3 unused, for 16 byte alignment only
//
//4
//5
//6
//7 unused

//use COL0,COL1 to index vectors of Vec6b

//spatial matrices are 6x6 (four 3x3 blocks)
// 0,  4,  8,  24, 28, 32
// 1,  5,  9,  25, 29, 33
// 2,  6, 10,  26, 30, 34
// 3,  7, 11,  27, 31, 35 unused, for 16 byte alignment only
//
//12, 16, 20,  36, 40, 44
//13, 17, 21,  37, 41, 45
//14, 18, 22,  38, 42, 46
//15, 19, 23,  39, 43, 47 unused

#define BLOCK0 0
#define BLOCK1 12
#define BLOCK2 24
#define BLOCK3 36

//SET FUNCTIONS

//inline void setVec6b(const Real val0, const Real val1, const Real val2, const Real val3, const Real val4, const Real val5, Vec6b Dest) {
#define setVec6b(val0,val1,val2,val3,val4,val5,Dest) do { \
	setVec3(val0,val1,val2,Dest); \
	setVec3(val3,val4,val5,Dest+COL1); \
} while(0)


//inline void setVec6b(const Real val, Vec6b Dest) {
#define setcVec6b(val,Dest) do { \
	setVec3(val,val,val,Dest); \
	setVec3(val,val,val,Dest+COL1); \
} while(0)

//COPY FUNCTIONS

//inline void copyVec6b(const Vec6b Source, Vec6b Dest) {
#define copyVec6b(Source,Dest) do { \
	copyVec3(Source,Dest); \
	copyVec3(Source+COL1,Dest+COL1); \
} while(0)

//inline void setMat6b(const Real val, Mat6b Dest) {
#define setcMat6b( val, Dest ) do { \
	setcMat3(val,Dest); \
	setcMat3(val,Dest+BLOCK1); \
	setcMat3(val,Dest+BLOCK2); \
	setcMat3(val,Dest+BLOCK3); \
} while(0)

//inline void copyMat6b(const Mat6b Source, Mat6b Dest) {
#define copyMat6b(Source,Dest) do { \
	copyMat3(Source,Dest); \
	copyMat3(Source+BLOCK1,Dest+BLOCK1); \
	copyMat3(Source+BLOCK2,Dest+BLOCK2); \
	copyMat3(Source+BLOCK3,Dest+BLOCK3); \
} while(0)

//inline void copyTMat6b(const Mat6b P, Mat6b PT) {
#define copyTMat6b(P,PT) do { \
	copyTMat3(P,PT); \
	copyTMat3(P+BLOCK1,PT+BLOCK2); \
	copyTMat3(P+BLOCK2,PT+BLOCK1); \
	copyTMat3(P+BLOCK3,PT+BLOCK3); \
} while(0)


//ARITHMETIC FUNCTIONS
//ADD

#if LINALG3_SIMD
//TODO, use for loops for auto-vectorization
#else

//inline void addVec6b(const Vec3 A, const Vec3 B, Vec3 C) {
#define addVec6b(A,B,C) do { \
	addVec3(A,B,C); \
	addVec3(A+COL1,B+COL1,C+COL1); \
} while(0)


//inline void addMat6b(const Vec3 A, const Vec3 B, Vec3 C) {
#define addMat6b(A,B,C) do { \
	addMat3(A,B,C); \
	addMat3(A+BLOCK1,B+BLOCK1,C+BLOCK1); \
	addMat3(A+BLOCK2,B+BLOCK2,C+BLOCK2); \
	addMat3(A+BLOCK3,B+BLOCK3,C+BLOCK3); \
} while(0)

#endif


//CROSS PRODUCT OPERAITONS
//cross product with motion vector (m)
//d = v x m
//[d0] = [skew(v0) 0	   ] [m0]
//[d1]   [skew(v1) skew(v0)] [m1]
//d0 = v0 x m0
//d1 = v1 x m0 + v0 x m1
//inline void crossVec6bMotion(const Vec6b v, const Vec6b m, Vec6b d) {
#define crossVec6bMotion(v,m,d) do { \
	crossVec3(v,m+COL1,d); /*temporary*/ \
	crossVec3(v+COL1,m,d+COL1); \
	addVec3(d,d+COL1,d+COL1); \
	crossVec3(v,m,d); \
} while(0)

//cross product with force vector (f)
//d = v x f
//[d0] = [skew(v0) skew(v1)] [f0]
//[d1]   [0		   skew(v0)] [f1]
//d0 = v0 x f0 + v1 x f1
//d1 =			 v0 x f1
//inline void crossVec6bForce(const Vec6b v, const Vec6b f, Vec6b d) {
#define crossVec6bForce(v,f,d) do { \
	crossVec3(v,f,d+COL1); /*temporary*/ \
	crossVec3(v+COL1,f+COL1,d); \
	addVec3(d,d+COL1,d); \
	crossVec3(v,f+COL1,d+COL1); \
} while(0)

//convert homogeneous transform to Plucker transform
//HT=[R|t]
//P=
//[R		 0]
//[skew(t)*R R]
inline void HTToPlucker(const HomogeneousTransform T, Mat6b P) {
//#define HTToPlucker(T,P) do { 
	skewVec3(T+COL3,P); /*temporary*/ \
	multMatMat3(P,T,P+BLOCK1); \
	copyMat3(T,P); \
	copyMat3(T,P+BLOCK3); \
	setcMat3(0,P+BLOCK2); \
} //while(0)

	
//convert inverse of homogeneous transform to Plucker transform
//HT=[R|t]
//P=
//[ R'			0 ]
//[-R'*skew3(t) R']
inline void invHTToPlucker(const HomogeneousTransform T, Mat6b P) {
//#define invHTToPlucker(T,P) do { 
	skewVec3(T+COL3,P); /*temporary*/ \
	multMatTMat3(T,P,P+BLOCK1); \
	mulcMat3(P+BLOCK1,-1,P+BLOCK1); \
	copyTMat3(T,P); \
	copyTMat3(T,P+BLOCK3); \
	setcMat3(0,P+BLOCK2); \
} //while(0)


//convert Plucker transform to homogeneous transform
//HT = [P0 | unskew(P1*P0')]
inline void PluckerToHT(const Mat6b P, HomogeneousTransform T) {
	Mat3 Tmp;
	copyTMat3(P,Tmp);
	multMatMat3(P+BLOCK1,Tmp,T); //temporary
	unskewVec3(T,T+COL3);
	copyMat3(P,T);
}


//macros for Plucker transform
//use Plucker block 2 = zeros to reduce computation

//multiply spatial vector by Plucker transform
//w=P*v
//[w0] = [P0 0 ] [v0]
//[w1]   [P1 P3] [v1]
//w0 = P0*v0
//w1 = P1*v0 + P3*v1
//inline void multPluckerVec6b(const Mat6b P, const Vec6b v, Vec6b w) {
#define multPluckerVec6b(P,v,w) do { \
	multMatVec3(P+BLOCK1,v,w); /*temporary*/ \
	multMatVec3(P+BLOCK3,v+COL1,w+COL1); \
	addVec3(w,w+COL1,w+COL1); \
	multMatVec3(P,v,w); \
} while(0)

//multiply spatial vector by transpose of Plucker transform
//w=P'*v
//[w0] = [P0' P1'] [v0]
//[w1]   [0   P3'] [v1]
//w0 = P0'*v0 + P1'*v1
//w1 = P3'*v1
//inline void multPluckerTVec6b(const Mat6b P, const Vec6b v, Vec6b w) {
#define multPluckerTVec6b(P,v,w) do { \
	multMatTVec3(P,v,w+COL1); /*temporary*/ \
	multMatTVec3(P+BLOCK1,v+COL1,w); \
	addVec3(w,w+COL1,w); \
	multMatTVec3(P+BLOCK3,v+COL1,w+COL1); \
} while(0)

//multiply spatial vector by Mat6b
//w = M*v
//[w0] = [M0 M2] [v0]
//[w1]   [M1 M3] [v1]
//w0 = M0*v0 + M2*v1
//w1 = M1*v0 + M3*v1
inline void multMatVec6b(const Mat6b M, const Vec6b v, Vec6b w) {
	Vec3 tmp;
	multMatVec3(M,v,tmp);
	multMatVec3(M+BLOCK2,v+COL1,w);
	addVec3(w,tmp,w);
	multMatVec3(M+BLOCK1,v,tmp);
	multMatVec3(M+BLOCK3,v+COL1,w+COL1);
	addVec3(w+COL1,tmp,w+COL1);
}


//DON'T USE THESE FUNCTIONS DIRECTLY
//R=P'*M

//[R0 R2] = [P0' P1'] [M0 M2]
//[R1 R3]   [0   P3'] [M1 M3]
//R0 = P0'*M0 + P1'*M1
//R1 = P3'*M1;
//R2 = P0'*M2 + P1'*M3
//R3 = P3'*M3
inline void multPluckerTMat6b(const Mat6b P, const Mat6b M, Mat6b R) {
	Mat3 Tmp;

	multMatTMat3(P,M,Tmp);
	multMatTMat3(P+BLOCK1,M+BLOCK1,R);
	addMat3(R,Tmp,R);
	multMatTMat3(P+BLOCK3,M+BLOCK1,R+BLOCK1);
	multMatTMat3(P,M+BLOCK2,Tmp);
	multMatTMat3(P+BLOCK1,M+BLOCK3,R+BLOCK2);
	addMat3(R+BLOCK2,Tmp,R+BLOCK2);
	multMatTMat3(P+BLOCK3,M+BLOCK3,R+BLOCK3);
}


//R=M*P

//[R0 R2] = [M0 M2] [P0 0 ]
//[R1 R3]   [M1 M3] [P1 P3]
//R0 = M0*P0 + M2*P1;
//R1 = M1*P0 + M3*P1;
//R2 = M2*P3;
//R3 = M3*P3;
inline void multMat6bPlucker(const Mat6b M, const Mat6b P, Mat6b R) {
	Mat3 Tmp;

	multMatMat3(M,P,Tmp); 
	multMatMat3(M+BLOCK2,P+BLOCK1,R); 
	addMat3(R,Tmp,R); 
	multMatMat3(M+BLOCK1,P,Tmp); 
	multMatMat3(M+BLOCK3,P+BLOCK1,R+BLOCK1); 
	addMat3(R+BLOCK1,Tmp,R+BLOCK1); 
	multMatMat3(M+BLOCK2,P+BLOCK3,R+BLOCK2); 
	multMatMat3(M+BLOCK3,P+BLOCK3,R+BLOCK3); 
}



//USE THIS FUNCTION DIRECTLY
//transform matrix using Plucker transform
//R = P*M*P'
inline void multPluckerTMat6bPlucker( const Mat6b P, const Mat6b M, Mat6b R) {
	Mat6b Tmp;

	multMat6bPlucker(M,P,Tmp);
	multPluckerTMat6b(P,Tmp,R);
}


//copy to/from array without extra elements for alignment
inline void copyVec6bToArray(const Vec6b Source, Real* Dest) {
	copyVec3ToArray(2,(const Vec3*) Source, Dest);
}

inline void copyArrayToVec6b(const Real* Source, Vec6b Dest) {
	copyArrayToVec3(Source,2,(Vec3*) Dest);
}

inline void copyMat6bColToVec6b(const int ci, const Mat6b M, Vec6b Dest) {
	//ci: column index
	int j = (ci/3)*BLOCK2+(ci%3)*SIZEVEC3; //index of first element of column
	const Real* p = M+j; //pointer to first element of column
	setVec6b(	p[0],
				p[1],
				p[2],
				p[0+BLOCK1],
				p[1+BLOCK1],
				p[2+BLOCK1], Dest);
}

//TODO, eliminate duplication?
inline void copyMat6bColToArray(const int ci, const Mat6b M, Real* Dest) {
	//ci: column index
	int j = (ci/3)*BLOCK2+(ci%3)*SIZEVEC3; //index of first element of column
	const Real* p = M+j; //pointer to first element of column
	Dest[0]=p[0];
	Dest[1]=p[1];
	Dest[2]=p[2];
	Dest[3]=p[0+BLOCK1];
	Dest[4]=p[1+BLOCK1];
	Dest[5]=p[2+BLOCK1];
}

inline void copyMat6bToArray(const Mat6b Source, Real* Dest) {
	for (int j=0; j<6; j++) //loop over columns
		copyMat6bColToArray(j,Source,Dest+S2I(0,j,6));
}


//PRINT FUNCTIONS
inline void printVec6b(const Vec6b V, int precision, int width) {
	Real V_[6];
	copyVec3ToArray(2,(const Vec3*) V, V_);
	printMatReal(6,1,V_,precision,width);
}

inline void printnVec6b(const int n, const Vec6b V[], int precision, int width) { 

	Real* V_;
	V_ = new Real[6*n];

	copyVec3ToArray(2*n,(const Vec3*) V, V_);
	printMatReal(6,n,V_,precision,width); 

	delete[] V_;
}

inline void printMat6b(const Mat6b M, int precision, int width) {
	Real M_[6*6];
	copyMat6bToArray(M,M_);
	printMatReal(6,6,M_,precision,width);
}

void toSpatialInertia(const Real m, const Vec3 c, const Mat3 I, Mat6b Is);
void fromSpatialInertia(const Mat6b Is, Real m, Vec3 c, Mat3 I);

void multPluckerTInertiaPlucker(const Mat6b P, const Mat6b I, Mat6b R); //faster alternative

#endif //_WMRSIM_SPATIAL_H_