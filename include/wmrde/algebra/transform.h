//utility for Homogeneous Transforms

#ifndef _WMRDE_TRANSFORM_H_
#define _WMRDE_TRANSFORM_H_

#include <iostream>
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
  inline HTransform inverse() const
  {
    return HTransform(R.transpose(), -R.transpose()*t);
  }

  inline HTransform composeWith(const HTransform& other) const
  {
    //return (*this) * other
    return HTransform( R*other.R, R*other.t + t);

    //TODO, is copy elision performed?
    //TODO, is this faster?
//    HTransform out;
//    out.R.noalias() = R*other.R;
//    out.t.noalias() = R*other.t + t;
//    return out;
  }

  inline HTransform composeInvWith(const HTransform& other) const
  {
    //return (*this.inverse()) * other
    return HTransform( R.transpose()*other.R, R.transpose()*(other.t - t));
  }

  inline Vec3 applyTo(const Vec3& vec) const
  {
    return R*vec + t;
  }
  inline Vec3 applyInvTo(const Vec3& vec) const
  {
    return R.transpose()*(vec - t);
  }

  //FOR DEBUGGING
  inline Eigen::Matrix<Real,3,4> to3x4() const
  {
    Eigen::Matrix<Real,3,4> out;
    out.block(0,0,3,3) = R;
    out.block(0,3,3,1) = t;
    return out;
  }

  inline Eigen::Matrix<Real,4,4> to4x4() const
  {
    Eigen::Matrix<Real,4,4> out;
    out.setZero();
    out.block(0,0,3,3) = R;
    out.block(0,3,3,1) = t;
    out(3,3) = 1.0;
    return out;
  }

  friend std::ostream &operator<<( std::ostream &output, const HTransform &T)
  {
    output << T.to3x4();
    return output;
  }

  inline bool isApprox(const HTransform& other, const Real prec = 0) const
  {
    return this->to3x4().isApprox(other.to3x4(), prec);
  }
};

} //namespace

#endif //_WMRDE_TRANSFORM_H_
