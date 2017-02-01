//utility for Homogeneous Transforms

#ifndef _WMRDE_TRANSFORM_H_
#define _WMRDE_TRANSFORM_H_

#include <wmrde/algebra/linalg3.h>

namespace wmrde
{

class HTransform //homogeneous transform
{
public:
  Mat3 R; //Rotation matrix
  Vec3 t; //translation vector

  HTransform() {};
  HTransform(const Mat3& Rotation, const Vec3& translation) :
    R(Rotation),
    t(translation)
  {}

  inline void setIdentity() { R.setIdentity(); t.setZero(); }
  inline void invert() //invert in place
  {
    R.transposeInPlace();
    t = -R*t; //here R is transpose of original R
  }

  inline HTransform compose(const HTransform& other)
  {
    //return (*this) * other
    return HTransform( R*other.R, R*other.t + t);
  }

  inline HTransform composeInv(const HTransform& other)
  {
    //return (*this.inverse()) * other
    return HTransform( R.transpose()*other.R, R.transpose()*(other.t - t));
  }

  inline Vec3 apply(const Vec3& vec)
  {
    return R*vec + t;
  }
  inline Vec3 applyInv(const Vec3& vec)
  {
    return R.transpose()*(vec - t);
  }

  //DEBUGGING
  inline Eigen::Matrix<Real,3,4> to3x4()
  {
    Eigen::Matrix<Real,3,4> out;
    out.block(0,0,3,3) = R;
    out.block(0,3,3,1) = t;
    return out;
  }

  inline Eigen::Matrix<Real,4,4> to4x4()
  {
    Eigen::Matrix<Real,4,4> out;
    out.setZero();
    out.block(0,0,3,3) = R;
    out.block(0,3,3,1) = t;
    out(3,3) = 1.0;
    return out;
  }
};

} //namespace

#endif //_WMRDE_TRANSFORM_H_
