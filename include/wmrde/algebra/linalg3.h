//linalg3.h
//linear algebra for vectors/matrices in R3
//Vec3 has an extra element for 16 byte alignment for vectorization
//vectorized functions can only be used on Vec3 & Mat3 arrays!
//do not use on sub arrays with no extra element
//Visual Studio vectorizer messages:
//http://msdn.microsoft.com/en-us/library/jj658585.aspx


#ifndef _WMRDE_LINALG3_H_
#define _WMRDE_LINALG3_H_

#include <cmath> //for sqrt

#define LINALG3_SIMD 0

//array sizes
#define SIZEVEC3 4
#define SIZEMAT3 12

#include <wmrde/algebra/matrix.h>

typedef Real Vec3[SIZEVEC3];
typedef Real Mat3[SIZEMAT3];

//can convert between macros and inline functions
//macro pros: guaranteed to inline, faster?
//inline function pros: type safe, lets compiler decide whether or not to inline, can overload (but don't)


//Vec3 indices:
//0
//1
//2
//3 unused, for 16 byte alignment only

//Mat3 indices: 
//(column-major order)
//0,4,8
//1,5,9
//2,6,10
//3,7,11 unused, for 16 byte alignment only

#define COL0 0
#define COL1 4
#define COL2 8

//don't vectorize set or copy functions, code 1300 no computation
//safe to use on sub arrays

//SET FUNCTIONS
#define setVec3(val0,val1,val2,Dest) do { \
	(Dest)[0] = (val0); \
	(Dest)[1] = (val1); \
	(Dest)[2] = (val2); \
} while(0)

#define setcVec3(val,Dest) do { \
	(Dest)[0] = (val); \
	(Dest)[1] = (val); \
	(Dest)[2] = (val); \
} while(0)

#define setMat3(val00,val10,val20,val01,val11,val21,val02,val12,val22,Dest) do { \
	setVec3(val00,val10,val20,Dest); \
	setVec3(val01,val11,val21,Dest+COL1); \
	setVec3(val02,val12,val22,Dest+COL2); \
} while(0)

#define setcMat3(val,Dest) \
	setMat3(val,val,val, val,val,val, val,val,val, Dest)

#define setDiagMat3(val00,val11,val22,Dest) \
	setMat3(val00,0,0, 0,val11,0, 0,0,val22, Dest)

#define setIdentityMat3(Dest) \
	setMat3(1.0,0,0, 0,1.0,0, 0,0,1.0, Dest)

//COPY FUNCTIONS
#define copyVec3(Source,Dest) do { \
	(Dest)[0] = (Source)[0]; \
	(Dest)[1] = (Source)[1]; \
	(Dest)[2] = (Source)[2]; \
} while(0)

#define copyMat3(Source,Dest) do { \
	copyVec3(Source,Dest); \
	copyVec3(Source+COL1,Dest+COL1); \
	copyVec3(Source+COL2,Dest+COL2); \
} while(0)

//TRANSPOSE FUNCTIONS
//pointer AT can't equal A
#define copyTMat3(A,AT) \
setMat3((A)[0],(A)[COL1],(A)[COL2], \
		(A)[1],(A)[1+COL1],(A)[1+COL2], \
		(A)[2],(A)[2+COL1],(A)[2+COL2], AT)

//copy a Mat3 row to Vec3
#define copyRowToVec3(ri,M,V) do { \
	(V)[0] = (M)[ri]; \
	(V)[1] = (M)[ri+COL1]; \
	(V)[2] = (M)[ri+COL2]; \
} while(0)

//copy a Vec3 to a Mat3 row
#define copyVec3ToRow(V,ri,M) do { \
	(M)[ri] = (V)[0]; \
	(M)[ri+COL1] = (V)[1]; \
	(M)[ri+COL2] = (V)[2]; \
} while(0)


//ARITHMETIC FUNCTIONS
//pointer C can't equal A or B
#if LINALG3_SIMD
//TODO, use for loops for auto-vectorization

#else
//Vec3
//ADD
#define addVec3(A,B,C) do { \
	(C)[0] = (A)[0] + (B)[0]; \
	(C)[1] = (A)[1] + (B)[1]; \
	(C)[2] = (A)[2] + (B)[2]; \
} while(0)

//add with a multiplier, m = -1 to subtract
#define addmVec3(A,m,B,C) do { \
	(C)[0] = (A)[0] + (m)*(B)[0]; \
	(C)[1] = (A)[1] + (m)*(B)[1]; \
	(C)[2] = (A)[2] + (m)*(B)[2]; \
} while(0)

//add a constant
#define addcVec3(A,b,C) do { \
	(C)[0] = (A)[0] + (b); \
	(C)[1] = (A)[1] + (b); \
	(C)[2] = (A)[2] + (b); \
} while(0)


//MULTIPLY
//inline void mulcVec3(const Vec3 A, const Real b, Vec3 C) {
#define mulcVec3(A,b,C) do { \
	(C)[0] = (A)[0]*(b); \
	(C)[1] = (A)[1]*(b); \
	(C)[2] = (A)[2]*(b); \
} while(0)

//Mat3
//ADD
//inline void addMat3(const Mat3 A, const Mat3 B, Mat3 C) {
#define addMat3(A,B,C) do { \
	addVec3(A,B,C); \
	addVec3(A+COL1,B+COL1,C+COL1); \
	addVec3(A+COL2,B+COL2,C+COL2); \
} while(0)

#define addmMat3(A,m,B,C) do { \
	addmVec3(A,m,B,C); \
	addmVec3(A+COL1,m,B+COL1,C+COL1); \
	addmVec3(A+COL2,m,B+COL2,C+COL2); \
} while(0)

#define addcMat3(A,b,C) do { \
	addcVec3(A,b,C); \
	addcVec3(A+COL1,b,C+COL1); \
	addcVec3(A+COL2,b,C+COL2); \
} while(0)


//MULTIPLY
#define mulcMat3(A,b,C) do { \
	mulcVec3(A,b,C); \
	mulcVec3(A+COL1,b,C+COL1); \
	mulcVec3(A+COL2,b,C+COL2); \
} while(0)

#endif


//VECTOR OPERATIONS
#define normVec3(A) sqrt((A)[0]*(A)[0] + (A)[1]*(A)[1] + (A)[2]*(A)[2])
//must use inline function to precompute normal
inline void normalizeVec3(Vec3 A) {
	Real n = normVec3(A);
	mulcVec3(A,1.0/n,A);
}
#define dotVec3(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])


//C = A x B
//pointer C can't equal A or B
#define crossVec3(A,B,C) do { \
	(C)[0] = -(A)[2]*(B)[1] + (A)[1]*(B)[2]; \
	(C)[1] =  (A)[2]*(B)[0] - (A)[0]*(B)[2]; \
	(C)[2] = -(A)[1]*(B)[0] + (A)[0]*(B)[1]; \
} while(0)

//A x B = skew(A)*B
//  0, -v2,  v1
// v2,   0, -v0
//-v1,  v0,   0 
#define skewVec3(V,M) do { \
	(M)[0]= 0;		(M)[0+COL1]=-(V)[2];	(M)[0+COL2]= (V)[1]; \
	(M)[1]= (V)[2];	(M)[1+COL1]= 0;			(M)[1+COL2]=-(V)[0]; \
	(M)[2]=-(V)[1];	(M)[2+COL1]= (V)[0];	(M)[2+COL2]= 0; \
} while(0)

#define unskewVec3(M,V) do { \
	(V)[0] = (M)[6]; \
	(V)[1] =-(M)[2]; \
	(V)[2] = (M)[1]; \
} while(0)



//MATRIX MULTIPLICATION
//c = A*b
//pointer c can't equal b
#if LINALG3_SIMD
//DEBUGGING, to vectorize

inline void multMatVec3(const Mat3 A, const Vec3 b, Vec3 c) {
	for (int i=0; i<NVEC3LOOP; i++)
		c[i] = A[i]*b[0] + (A+COL1)[i]*b[1] + (A+COL2)[i]*b[2]; //not vectorized code 1200
}

#else

//inline void multMatVec3(const Mat3 A, const Vec3 b, Vec3 c) {
#define multMatVec3(A,b,c) do { \
	(c)[0] = (A)[0]*(b)[0] + (A)[0+COL1]*(b)[1] + (A)[0+COL2]*(b)[2]; \
	(c)[1] = (A)[1]*(b)[0] + (A)[1+COL1]*(b)[1] + (A)[1+COL2]*(b)[2]; \
	(c)[2] = (A)[2]*(b)[0] + (A)[2+COL1]*(b)[1] + (A)[2+COL2]*(b)[2]; \
} while(0)

#endif



//c = (A^T)*b
//possible to vectorize?
//inline void multMatTVec3(const Mat3 A, const Vec3 b, Vec3 c) { 
#define multMatTVec3(A,b,c) do { \
	(c)[0] = dotVec3(A,b); \
	(c)[1] = dotVec3(A+COL1,b); \
	(c)[2] = dotVec3(A+COL2,b); \
} while(0)


//C = A*B
//inline void multMatMat3(const Mat3 A, const Mat3 B, Mat3 C) { 
#define multMatMat3(A,B,C) do { \
	multMatVec3(A,B,C); \
	multMatVec3(A,B+COL1,C+COL1); \
	multMatVec3(A,B+COL2,C+COL2); \
} while(0)

//C = (A^T)*B
//inline void multMatTMat3(const Mat3 A, const Mat3 B, Mat3 C) {
#define multMatTMat3(A,B,C) do { \
	multMatTVec3(A,B,C); \
	multMatTVec3(A,B+COL1,C+COL1); \
	multMatTVec3(A,B+COL2,C+COL2); \
} while(0)



//copy to/from array without extra elements for alignment
inline void copyVec3ToArray(const int n, const Vec3 Source[], Real* Dest) {
	for (int i=0; i<n; i++)
		copyVec3(Source[i],Dest+(i*3));
}

inline void copyArrayToVec3(const Real* Source, const int n, Vec3 Dest[]) {
	for (int i=0; i<n; i++)
		copyVec3(Source+(i*3),Dest[i]);
}

//PRINT FUNCTIONS
inline void printVec3(const Vec3 V, int precision, int width) { 
	printMatReal(3,1,V,precision,width); 
}

inline void printnVec3(const int n, const Vec3 V[], int precision, int width) {
	//dynamic memory allocation required
	Real* V_;
	V_ = new Real[n*3];

	copyVec3ToArray(n,V,V_);
	printMatReal(3,n,V_,precision,width);

	delete[] V_;
}

inline void printMat3(const Mat3 M, int precision, int width) { 
	Real M_[3*3];
	copyVec3ToArray(3,(const Vec3*) M,M_);
	printMatReal(3,3,M_,precision,width);
}

#endif //_WMRDE_LINALG3_H_
