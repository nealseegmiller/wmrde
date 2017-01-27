//linear algebra functions for vectors/matrices in R3

#ifndef _WMRDE_LINALG3_H_
#define _WMRDE_LINALG3_H_

#include <cmath> //for sqrt
#include <sstream>
#include <iomanip> //for setprecision

#include <wmrde/common.h>

namespace wmrde
{

#define VEC3_SIZE 4
#define MAT3_SIZE 12

typedef Real Vec3[VEC3_SIZE];
typedef Real Mat3[MAT3_SIZE]; //contains three Vec3 arrays

//Vec3 has an extra element for 16 byte alignment for SIMD vectorization.
//Vectorized functions can only be used on Vec3 & Mat3 arrays!
//Do not use on sub arrays with no extra element.

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

//TODO, undefine these at the end of the header file
#define MAT3_COL0 0 //first index in column 0 of Mat3
#define MAT3_COL1 4 //first index in column 1 of Mat3
#define MAT3_COL2 8 //first index in column 2 of Mat3


//SET FUNCTIONS
inline void setVec3(
    const Real val, Vec3 dst)
{
  dst[0] = val;
  dst[1] = val;
  dst[2] = val;
}

inline void setVec3(
    const Real val0, const Real val1, const Real val2, Vec3 dst)
{
	dst[0] = val0;
	dst[1] = val1;
	dst[2] = val2;
}

//specify elements in column-major order
inline void setMat3(
    const Real val00, const Real val10, const Real val20, //column 0
    const Real val01, const Real val11, const Real val21, //column 1
    const Real val02, const Real val12, const Real val22, //column 2
    Mat3 dst)
{
  dst[0] = val00;
  dst[1] = val10;
  dst[2] = val20;
//  dst[3] = 0.0; //unused
  dst[4] = val01;
  dst[5] = val11;
  dst[6] = val21;
//  dst[7] = 0.0; //unused
  dst[8] = val02;
  dst[9] = val12;
  dst[10] = val22;
//  dst[11] = 0.0; //unused
}

inline void setMat3(
    const Real val, Mat3 dst)
{
	setMat3(val,val,val, val,val,val, val,val,val, dst);
//  for (size_t i=0; i < MAT3_SIZE; i++) { dst[i] = val; } //TODO, faster?
}

//set to diagonal matrix
inline void setMat3Diag(
    const Real val00, const Real val11, const Real val22,
    Mat3 dst)
{
  setMat3(val00,0,0, 0,val11,0, 0,0,val22, dst);
}

//set to identity matrix
inline void setMat3Identity(Mat3 dst)
{
	setMat3(1.0,0,0, 0,1.0,0, 0,0,1.0, dst);
}

//COPY FUNCTIONS
inline void copyVec3(const Vec3 src, Vec3 dst)
{
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
}

//copy 3 elements from one array to another
inline void copy3(const Real* src, Real* dst)
{
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];
}

inline void copyMat3(const Mat3 src, Mat3 dst)
{
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];

  dst[4] = src[4];
  dst[5] = src[5];
  dst[6] = src[6];

  dst[8] = src[8];
  dst[9] = src[9];
  dst[10] = src[10];

//  for (size_t i=0; i < MAT3_SIZE; i++) { dst[i] = src[i]; } //TODO, faster?
}

//TRANSPOSE FUNCTIONS
//src and dst can't be the same (can't transpose in place)
inline void copyTMat3(const Mat3 src, Mat3 dst)
{
  dst[0] = src[0];
  dst[1] = src[4];
  dst[2] = src[8];

  dst[4] = src[1];
  dst[5] = src[5];
  dst[6] = src[9];

  dst[8] = src[2];
  dst[9] = src[6];
  dst[10] = src[10];
}

//copy a Mat3 row to Vec3
inline void copyMat3Row(const Mat3 src, const int row_idx, Vec3 dst)
{
	dst[0] = src[row_idx + MAT3_COL0];
	dst[1] = src[row_idx + MAT3_COL1];
	dst[2] = src[row_idx + MAT3_COL2];
}

//copy to/from Mat3 and an array without alignment elements
//useful for converting to/from Eigen matrices
inline void copyMat3ToArray(const Mat3 src, Real* dst)
{
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];

  dst[3] = src[4];
  dst[4] = src[5];
  dst[5] = src[6];

  dst[6] = src[8];
  dst[7] = src[9];
  dst[8] = src[10];
}

inline void copyArrayToMat3(const Real* src, Mat3 dst)
{
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];

  dst[4] = src[3];
  dst[5] = src[4];
  dst[6] = src[5];

  dst[8] = src[6];
  dst[9] = src[7];
  dst[10] = src[8];
}


inline bool Vec3Equal(const Mat3 mat_a, const Mat3 mat_b)
{
  Real tol = 10.0*std::numeric_limits<Real>::epsilon();
  bool equal = true;
  equal &= std::abs(mat_a[0] - mat_b[0]) < tol;
  equal &= std::abs(mat_a[1] - mat_b[1]) < tol;
  equal &= std::abs(mat_a[2] - mat_b[2]) < tol;

  return equal;
}

inline bool Mat3Equal(const Mat3 mat_a, const Mat3 mat_b)
{
  Real tol = 10.0*std::numeric_limits<Real>::epsilon();
  bool equal = true;
  equal &= std::abs(mat_a[0] - mat_b[0]) < tol;
  equal &= std::abs(mat_a[1] - mat_b[1]) < tol;
  equal &= std::abs(mat_a[2] - mat_b[2]) < tol;

  equal &= std::abs(mat_a[4] - mat_b[4]) < tol;
  equal &= std::abs(mat_a[5] - mat_b[5]) < tol;
  equal &= std::abs(mat_a[6] - mat_b[6]) < tol;

  equal &= std::abs(mat_a[8] - mat_b[8]) < tol;
  equal &= std::abs(mat_a[9] - mat_b[9]) < tol;
  equal &= std::abs(mat_a[10] - mat_b[10]) < tol;

  return equal;
}

//ARITHMETIC FUNCTIONS
#define LINALG3_SIMD 0
#if LINALG3_SIMD
//TODO, use for loops for auto-vectorization

#else

//VEC3 FUNCTIONS

//add vec_b to vec_a
inline void addVec3(const Vec3 vec_a, const Vec3 vec_b, Vec3 out)
{
	out[0] = vec_a[0] + vec_b[0];
	out[1] = vec_a[1] + vec_b[1];
	out[2] = vec_a[2] + vec_b[2];
}

//add m*vec_b to vec_a
inline void addmVec3(const Vec3 vec_a, const Real m, const Vec3 vec_b, Vec3 out)
{
  out[0] = vec_a[0] + m*vec_b[0];
  out[1] = vec_a[1] + m*vec_b[1];
  out[2] = vec_a[2] + m*vec_b[2];
}

//subtract vec_b from vec_a
//TODO, could instead use addmVec3(vec_a, -1.0, vec_b, out)
inline void subVec3(const Vec3 vec_a, const Vec3 vec_b, Vec3 out)
{
  out[0] = vec_a[0] - vec_b[0];
  out[1] = vec_a[1] - vec_b[1];
  out[2] = vec_a[2] - vec_b[2];
}

//multiply a vector by a coefficient
inline void mulcVec3(const Vec3 vec, const Real coeff, Vec3 out) {
	out[0] = coeff*vec[0];
	out[1] = coeff*vec[1];
	out[2] = coeff*vec[2];
}

//MAT3 FUNCTIONS

//add mat_b to mat_a
inline void addMat3(const Mat3 mat_a, const Mat3 mat_b, Mat3 out)
{
	out[0] = mat_a[0] + mat_b[0];
	out[1] = mat_a[1] + mat_b[1];
	out[2] = mat_a[2] + mat_b[2];

  out[4] = mat_a[4] + mat_b[4];
  out[5] = mat_a[5] + mat_b[5];
  out[6] = mat_a[6] + mat_b[6];

  out[8] = mat_a[8] + mat_b[8];
  out[9] = mat_a[9] + mat_b[9];
  out[10] = mat_a[10] + mat_b[10];
}

//add m*mat_b to mat_a
//to subtract use addmMat3(mat_a, -1.0, mat_b, out)
inline void addmMat3(const Mat3 mat_a, const Real m, const Mat3 mat_b, Mat3 out)
{
  out[0] = mat_a[0] + m*mat_b[0];
  out[1] = mat_a[1] + m*mat_b[1];
  out[2] = mat_a[2] + m*mat_b[2];

  out[4] = mat_a[4] + m*mat_b[4];
  out[5] = mat_a[5] + m*mat_b[5];
  out[6] = mat_a[6] + m*mat_b[6];

  out[8] = mat_a[8] + m*mat_b[8];
  out[9] = mat_a[9] + m*mat_b[9];
  out[10] = mat_a[10] + m*mat_b[10];
}

//multiply a matrix by a coefficient
inline void mulcMat3(const Mat3 mat, const Real coeff, Mat3 out) {
  out[0] = coeff*mat[0];
  out[1] = coeff*mat[1];
  out[2] = coeff*mat[2];

  out[4] = coeff*mat[4];
  out[5] = coeff*mat[5];
  out[6] = coeff*mat[6];

  out[8] = coeff*mat[8];
  out[9] = coeff*mat[9];
  out[10] = coeff*mat[10];
}

#endif


//VECTOR OPERATIONS
//return the norm of a vector
inline Real normVec3(const Vec3 vec)
{
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}
//normalize a vector
inline void normalizeVec3(Vec3 vec)
{
	mulcVec3(vec, 1.0/normVec3(vec), vec);
}
//compute the dot product of 2 vectors
inline Real dotVec3(const Vec3 vec_a, const Vec3 vec_b)
{
  return vec_a[0]*vec_b[0] + vec_a[1]*vec_b[1] + vec_a[2]*vec_b[2];
}

//compute the cross product of vec_a and vec_b
inline void crossVec3(const Vec3 vec_a, const Vec3 vec_b, Vec3 out)
{
	setVec3(
	    -vec_a[2]*vec_b[1] + vec_a[1]*vec_b[2],
	     vec_a[2]*vec_b[0] - vec_a[0]*vec_b[2],
	    -vec_a[1]*vec_b[0] + vec_a[0]*vec_b[1], out);
}

//compute skew-symmetric matrix from vector for computing cross products
//vec_a x vec_b = skew(vec_a)*vec_b
//skew(vec_a) = [
//  0, -v2,  v1
// v2,   0, -v0
//-v1,  v0,   0]
inline void skewVec3(const Vec3 vec, Mat3 out)
{
  setMat3(0, vec[2], -vec[1],
      -vec[2], 0, vec[0],
       vec[1], -vec[0], 0, out);
}

inline void unskewVec3(const Mat3 mat, Vec3 out)
{
	out[0] = mat[6];
	out[1] =-mat[2];
	out[2] = mat[1];
}

//MATRIX MULTIPLICATION
//TODO, possible to vectorize?

//out = mat*vec
//out can't be the same as vec (can't multiply in place)
//TODO, use setVec3 to allow multiply in place?
inline void multMatVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
{
  out[0] = mat[0]*vec[0] + mat[4]*vec[1] + mat[8]*vec[2];
  out[1] = mat[1]*vec[0] + mat[5]*vec[1] + mat[9]*vec[2];
	out[2] = mat[2]*vec[0] + mat[6]*vec[1] + mat[10]*vec[2];
}

//out = (mat^T)*vec
inline void multMatTVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
{
	out[0] = mat[0]*vec[0] + mat[1]*vec[1] + mat[2]*vec[2];
	out[1] = mat[4]*vec[0] + mat[5]*vec[1] + mat[6]*vec[2];
	out[2] = mat[8]*vec[0] + mat[9]*vec[1] + mat[10]*vec[2];
}

//out = mat_a * mat_b
inline void multMatMat3(const Mat3 mat_a, const Mat3 mat_b, Mat3 out)
{
	multMatVec3(mat_a, &mat_b[MAT3_COL0], &out[MAT3_COL0]);
	multMatVec3(mat_a, &mat_b[MAT3_COL1], &out[MAT3_COL1]);
	multMatVec3(mat_a, &mat_b[MAT3_COL2], &out[MAT3_COL2]);
}

//out = (mat_a^T) * mat_b
inline void multMatTMat3(const Mat3 mat_a, const Mat3 mat_b, Mat3 out)
{
  multMatTVec3(mat_a, &mat_b[MAT3_COL0], &out[MAT3_COL0]);
  multMatTVec3(mat_a, &mat_b[MAT3_COL1], &out[MAT3_COL1]);
  multMatTVec3(mat_a, &mat_b[MAT3_COL2], &out[MAT3_COL2]);
}

//PRINT FUNCTIONS
//TODO, use wmrde::Matrix class for smarter formatting
inline std::string Vec3ToString(const Vec3 vec, const int precision = 6)
{
  std::stringstream ss;
  ss << std::setprecision(precision);
  ss << vec[0] << vec[1] << vec[2] << "\n";
  return ss.str();
}

inline std::string Mat3ToString(const Mat3 mat, const int precision = 6)
{
  std::stringstream ss;
  ss << std::setprecision(precision);
  ss << mat[0] << " " << mat[4] << " " << mat[8] << "\n"
     << mat[1] << " " << mat[5] << " " << mat[9] << "\n"
     << mat[2] << " " << mat[6] << " " << mat[10] << "\n";

  return ss.str();
}

} //namespace

#endif //_WMRDE_LINALG3_H_
