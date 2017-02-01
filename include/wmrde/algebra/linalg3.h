//linear algebra functions for vectors/matrices in R3

#ifndef _WMRDE_LINALG3_H_
#define _WMRDE_LINALG3_H_

#include <cmath> //for sqrt
#include <sstream>
#include <iomanip> //for setprecision
#include <limits>

#include <wmrde/common.h>

namespace wmrde
{

#define LINALG3_SIMD 0
//TODO, why are the vectorized implementations so much slower?

//#if LINALG3_SIMD
#if 1

#define VEC3_SIZE 4
#define MAT3_SIZE 12

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

enum Mat3Indices {
  MAT3_0=0,
  MAT3_1=1,
  MAT3_2=2,
  MAT3_3=4,
  MAT3_4=5,
  MAT3_5=6,
  MAT3_6=8,
  MAT3_7=9,
  MAT3_8=10
};

enum Mat3ColumnIndices {
  MAT3_COL0=0,
  MAT3_COL1=4,
  MAT3_COL2=8
};

#else

#define VEC3_SIZE 3
#define MAT3_SIZE 9

enum Mat3Indices {
  MAT3_0=0,
  MAT3_1=1,
  MAT3_2=2,
  MAT3_3=3,
  MAT3_4=4,
  MAT3_5=5,
  MAT3_6=6,
  MAT3_7=7,
  MAT3_8=8
};

enum Mat3ColumnIndices {
  MAT3_COL0=0,
  MAT3_COL1=3,
  MAT3_COL2=6
};

#endif

typedef Real Vec3[VEC3_SIZE];
typedef Real Mat3[MAT3_SIZE]; //contains three Vec3 arrays

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
  dst[MAT3_0] = val00;
  dst[MAT3_1] = val10;
  dst[MAT3_2] = val20;

  dst[MAT3_3] = val01;
  dst[MAT3_4] = val11;
  dst[MAT3_5] = val21;

  dst[MAT3_6] = val02;
  dst[MAT3_7] = val12;
  dst[MAT3_8] = val22;
}

inline void setMat3(
    const Real val, Mat3 dst)
{
  //set a continuous block of memory
  for (size_t i=0; i<MAT3_SIZE; i++) { dst[i] = val; }
}

//set to diagonal matrix
inline void setMat3Diagonal(
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

//DEBUGGING
inline void copy3n(const int n, const Real* src, Real* dst)
{
  for (int i = 0; i < n; i++)
  {
    copy3(src+(i*VEC3_SIZE), dst+(i*3));
  }
}

//TODO, this vectorizes but why is it slower?
//inline void copyMat3(const Mat3 src, Mat3 dst)
//{
//  //copy a contiguous block of memory
//  for (size_t i=0; i<MAT3_SIZE; i++) { dst[i] = src[i]; }
//}

inline void copyMat3(const Mat3 src, Mat3 dst)
{
  dst[MAT3_0] = src[MAT3_0];
  dst[MAT3_1] = src[MAT3_1];
  dst[MAT3_2] = src[MAT3_2];
  dst[MAT3_3] = src[MAT3_3];
  dst[MAT3_4] = src[MAT3_4];
  dst[MAT3_5] = src[MAT3_5];
  dst[MAT3_6] = src[MAT3_6];
  dst[MAT3_7] = src[MAT3_7];
  dst[MAT3_8] = src[MAT3_8];
}

#define copyMat3macro(src,dst) do { \
  dst[0] = src[0];\
  dst[1] = src[1];\
  dst[2] = src[2];\
  dst[4] = src[4];\
  dst[5] = src[5];\
  dst[6] = src[6];\
  dst[8] = src[8];\
  dst[9] = src[9];\
  dst[10]= src[10];\
} while(0)

//TRANSPOSE FUNCTIONS
//src and dst can't be the same (can't transpose in place)
inline void copyTMat3(const Mat3 src, Mat3 dst)
{
  //row 0 to col 0
  dst[MAT3_0] = src[MAT3_0];
  dst[MAT3_1] = src[MAT3_3];
  dst[MAT3_2] = src[MAT3_6];
  //row 1 to col 1
  dst[MAT3_3] = src[MAT3_1];
  dst[MAT3_4] = src[MAT3_4];
  dst[MAT3_5] = src[MAT3_7];
  //row 2 to col 2
  dst[MAT3_6] = src[MAT3_2];
  dst[MAT3_7] = src[MAT3_5];
  dst[MAT3_8] = src[MAT3_8];
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
  dst[0] = src[MAT3_0];
  dst[1] = src[MAT3_1];
  dst[2] = src[MAT3_2];

  dst[3] = src[MAT3_3];
  dst[4] = src[MAT3_4];
  dst[5] = src[MAT3_5];

  dst[6] = src[MAT3_6];
  dst[7] = src[MAT3_7];
  dst[8] = src[MAT3_8];
}

inline void copyArrayToMat3(const Real* src, Mat3 dst)
{
  dst[MAT3_0] = src[0];
  dst[MAT3_1] = src[1];
  dst[MAT3_2] = src[2];

  dst[MAT3_3] = src[3];
  dst[MAT3_4] = src[4];
  dst[MAT3_5] = src[5];

  dst[MAT3_6] = src[6];
  dst[MAT3_7] = src[7];
  dst[MAT3_8] = src[8];
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
  bool equal = true;
  equal &= Vec3Equal(&mat_a[MAT3_COL0], &mat_b[MAT3_COL0]);
  equal &= Vec3Equal(&mat_a[MAT3_COL1], &mat_b[MAT3_COL1]);
  equal &= Vec3Equal(&mat_a[MAT3_COL2], &mat_b[MAT3_COL2]);

  return equal;
}

//VEC3 ARITHMETIC FUNCTIONS
#if LINALG3_SIMD

//add vec_b to vec_a
inline void addVec3(const Vec3 vec_a, const Vec3 vec_b, Vec3 out)
{
  //This loop is unrolled but not vectorized
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] = vec_a[i] + vec_b[i]; }
}

//add m*vec_b to vec_a
inline void addmVec3(const Vec3 vec_a, const Real m, const Vec3 vec_b, Vec3 out)
{
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] = vec_a[i] + m*vec_b[i]; }
}

//subtract vec_b from vec_a
inline void subVec3(const Vec3 vec_a, const Vec3 vec_b, Vec3 out)
{
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] = vec_a[i] - vec_b[i]; }
}

//multiply vector by a scalar
inline void mulcVec3(const Vec3 vec, const Real c, Vec3 out) {
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] = c*vec[i]; }
}

#else

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

//multiply vector by a scalar
inline void mulcVec3(const Vec3 vec, const Real c, Vec3 out) {
  out[0] = c*vec[0];
  out[1] = c*vec[1];
  out[2] = c*vec[2];
}

#endif


//MAT3 ARITHMETIC FUNCTIONS

#if LINALG3_SIMD
//add mat_b to mat_a
inline void addMat3(const Mat3 mat_a, const Mat3 mat_b, Mat3 out)
{
  //This loop is vectorized
  for (size_t i=0; i<MAT3_SIZE; i++) { out[i] = mat_a[i] + mat_b[i]; }
}

//add m*mat_b to mat_a
//to subtract use addmMat3(mat_a, -1.0, mat_b, out)
inline void addmMat3(const Mat3 mat_a, const Real m, const Mat3 mat_b, Mat3 out)
{
  for (size_t i=0; i<MAT3_SIZE; i++) { out[i] = mat_a[i] + m*mat_b[i]; }
}

//multiply matrix by a scalar constant
inline void mulcMat3(const Mat3 mat, const Real c, Mat3 out)
{
  for (size_t i=0; i<MAT3_SIZE; i++) { out[i] = c*mat[i]; }
}

#else

//add mat_b to mat_a
inline void addMat3(const Mat3 mat_a, const Mat3 mat_b, Mat3 out)
{
	out[MAT3_0] = mat_a[MAT3_0] + mat_b[MAT3_0];
	out[MAT3_1] = mat_a[MAT3_1] + mat_b[MAT3_1];
	out[MAT3_2] = mat_a[MAT3_2] + mat_b[MAT3_2];

  out[MAT3_3] = mat_a[MAT3_3] + mat_b[MAT3_3];
  out[MAT3_4] = mat_a[MAT3_4] + mat_b[MAT3_4];
  out[MAT3_5] = mat_a[MAT3_5] + mat_b[MAT3_5];

  out[MAT3_6] = mat_a[MAT3_6] + mat_b[MAT3_6];
  out[MAT3_7] = mat_a[MAT3_7] + mat_b[MAT3_7];
  out[MAT3_8] = mat_a[MAT3_8] + mat_b[MAT3_8];
}

//add m*mat_b to mat_a
//to subtract use addmMat3(mat_a, -1.0, mat_b, out)
inline void addmMat3(const Mat3 mat_a, const Real m, const Mat3 mat_b, Mat3 out)
{
  out[MAT3_0] = mat_a[MAT3_0] + m*mat_b[MAT3_0];
  out[MAT3_1] = mat_a[MAT3_1] + m*mat_b[MAT3_1];
  out[MAT3_2] = mat_a[MAT3_2] + m*mat_b[MAT3_2];

  out[MAT3_3] = mat_a[MAT3_3] + m*mat_b[MAT3_3];
  out[MAT3_4] = mat_a[MAT3_4] + m*mat_b[MAT3_4];
  out[MAT3_5] = mat_a[MAT3_5] + m*mat_b[MAT3_5];

  out[MAT3_6] = mat_a[MAT3_6] + m*mat_b[MAT3_6];
  out[MAT3_7] = mat_a[MAT3_7] + m*mat_b[MAT3_7];
  out[MAT3_8] = mat_a[MAT3_8] + m*mat_b[MAT3_8];
}

//multiply matrix by a scalar constant
inline void mulcMat3(const Mat3 mat, const Real c, Mat3 out) {
  out[MAT3_0] = c*mat[MAT3_0];
  out[MAT3_1] = c*mat[MAT3_1];
  out[MAT3_2] = c*mat[MAT3_2];

  out[MAT3_3] = c*mat[MAT3_3];
  out[MAT3_4] = c*mat[MAT3_4];
  out[MAT3_5] = c*mat[MAT3_5];

  out[MAT3_6] = c*mat[MAT3_6];
  out[MAT3_7] = c*mat[MAT3_7];
  out[MAT3_8] = c*mat[MAT3_8];
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
#if LINALG3_SIMD

//These loops are vectorized but this method is slower!
inline void multMatVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
{
  setVec3(0.0, out);
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] += mat[MAT3_COL0+i]*vec[0]; }
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] += mat[MAT3_COL1+i]*vec[1]; }
  for (size_t i=0; i<VEC3_SIZE; i++) { out[i] += mat[MAT3_COL2+i]*vec[2]; }
}

//These loops are unrolled but not vectorized
inline void multMatTVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
{
  setVec3(0.0, out);
  for (size_t i=0; i<VEC3_SIZE; i++) { out[0] += mat[MAT3_COL0+i]*vec[i]; }
  for (size_t i=0; i<VEC3_SIZE; i++) { out[1] += mat[MAT3_COL1+i]*vec[i]; }
  for (size_t i=0; i<VEC3_SIZE; i++) { out[2] += mat[MAT3_COL2+i]*vec[i]; }
}

#else

//out = mat*vec
//out can't be the same as vec (can't multiply in place)
//TODO, use setVec3 to allow multiply in place?
//inline void multMatVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
//{
#define multMatVec3macro(mat, vec, out) do { \
  out[0] = mat[0]*vec[0] + mat[4]*vec[1] + mat[8]*vec[2]; \
  out[1] = mat[1]*vec[0] + mat[5]*vec[1] + mat[9]*vec[2]; \
	out[2] = mat[2]*vec[0] + mat[6]*vec[1] + mat[10]*vec[2]; \
} while(0)

//out = mat*vec
//out can't be the same as vec (can't multiply in place)
//TODO, use setVec3 to allow multiply in place?
inline void multMatVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
{
  out[0] = mat[MAT3_0]*vec[0] + mat[MAT3_3]*vec[1] + mat[MAT3_6]*vec[2];
  out[1] = mat[MAT3_1]*vec[0] + mat[MAT3_4]*vec[1] + mat[MAT3_7]*vec[2];
  out[2] = mat[MAT3_2]*vec[0] + mat[MAT3_5]*vec[1] + mat[MAT3_8]*vec[2];
}

//out = (mat^T)*vec
inline void multMatTVec3(const Mat3 mat, const Vec3 vec, Vec3 out)
{
	out[0] = mat[MAT3_0]*vec[0] + mat[MAT3_1]*vec[1] + mat[MAT3_2]*vec[2];
	out[1] = mat[MAT3_3]*vec[0] + mat[MAT3_4]*vec[1] + mat[MAT3_5]*vec[2];
	out[2] = mat[MAT3_6]*vec[0] + mat[MAT3_7]*vec[1] + mat[MAT3_8]*vec[2];
}

#endif

//out = mat_a * mat_b
inline void multMatMat3(const Mat3 mat_a, const Mat3 mat_b, Mat3 out)
{
	multMatVec3(mat_a, &mat_b[MAT3_COL0], &out[MAT3_COL0]);
	multMatVec3(mat_a, &mat_b[MAT3_COL1], &out[MAT3_COL1]);
	multMatVec3(mat_a, &mat_b[MAT3_COL2], &out[MAT3_COL2]);
}

#define multMatMat3macro(mat_a, mat_b, out) do { \
  multMatVec3(mat_a, &mat_b[0], &out[0]); \
  multMatVec3(mat_a, &mat_b[4], &out[4]); \
  multMatVec3(mat_a, &mat_b[8], &out[8]); \
} while(0)

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
  ss << vec[0] << "\n" << vec[1] << "\n" << vec[2] << "\n";
  return ss.str();
}

inline std::string Mat3ToString(const Mat3 mat, const int precision = 6)
{
  std::stringstream ss;
  ss << std::setprecision(precision);
  ss << mat[MAT3_0] << " " << mat[MAT3_3] << " " << mat[MAT3_6] << "\n"
     << mat[MAT3_1] << " " << mat[MAT3_4] << " " << mat[MAT3_7] << "\n"
     << mat[MAT3_2] << " " << mat[MAT3_5] << " " << mat[MAT3_8] << "\n";

  return ss.str();
}

} //namespace

#endif //_WMRDE_LINALG3_H_
