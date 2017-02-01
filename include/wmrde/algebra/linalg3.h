//utility for linear algebra in R3

#ifndef _WMRDE_LINALG3_H_
#define _WMRDE_LINALG3_H_

#include <wmrde/common.h>
#include <Eigen/Dense>

namespace wmrde
{

typedef Eigen::Matrix<Real,3,3> Mat3;
typedef Eigen::Matrix<Real,3,1> Vec3;

//compute skew-symmetric matrix from vector for computing cross products
//a.cross(b) = skew(a)*b
inline Mat3 skew(const Vec3 vec)
{
  Mat3 out;
  out << 0,     -vec(2), vec(1),
         vec(2), 0,     -vec(0),
        -vec(1), vec(0), 0;
  return out;
}

inline Vec3 unskew(const Mat3 mat)
{
  Vec3 out;
  out << mat(2,1), -mat(2,0), mat(1,0);
  return out;
}

} //namespace

#endif //_WMRDE_LINALG3_H_


