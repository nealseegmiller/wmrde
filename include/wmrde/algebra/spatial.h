//utility for spatial vector algebra for rigid body dynamics
//Reference: Roy Featherstone, Rigid Body Dynamics Algorithms
//http://royfeatherstone.org/spatial/index.html#spatial-software

#ifndef _WMRDE_SPATIAL_H_
#define _WMRDE_SPATIAL_H_

#include <wmrde/algebra/transform.h>

namespace wmrde
{

typedef Eigen::Matrix<Real,6,6> Mat6;
typedef Eigen::Matrix<Real,6,1> Vec6;


//CROSS PRODUCT OPERAITONS
//cross product with motion vector (m)
//d = v x m
//[d0] = [skew(v0) 0       ] [m0]
//[d1]   [skew(v1) skew(v0)] [m1]
//d0 = v0 x m0
//d1 = v1 x m0 + v0 x m1
inline Vec6 crossVecMotion(const Vec6& v, const Vec6& m)
{
  Vec6 out;
  out.topRows<3>() = v.topRows<3>().cross(m.topRows<3>());
  out.bottomRows<3>() = v.bottomRows<3>().cross(m.topRows<3>()) +
      v.topRows<3>().cross(m.bottomRows<3>());
  return out;
}

//cross product with force vector (f)
//d = v x f
//[d0] = [skew(v0) skew(v1)] [f0]
//[d1]   [0        skew(v0)] [f1]
//d0 = v0 x f0 + v1 x f1
//d1 = v0 x f1
inline Vec6 crossVecForce(const Vec6& v, const Vec6& f)
{
  Vec6 out;
  out.topRows<3>() = v.topRows<3>().cross(f.topRows<3>()) +
      v.bottomRows<3>().cross(f.bottomRows<3>());
  out.bottomRows<3>() = v.topRows<3>().cross(f.bottomRows<3>());
  return out;
}

//convert homogeneous transform to Plucker transform
//HT=[R|t]
//P=
//[R         0]
//[skew(t)*R R]
inline Mat6 HTToPlucker(const HTransform& HT)
{
  Mat6 out;
  out.topLeftCorner<3,3>() = HT.R;
  out.bottomLeftCorner<3,3>().noalias() = skew(HT.t)*HT.R;
  out.topRightCorner<3,3>().setZero();
  out.bottomRightCorner<3,3>() = HT.R;
  return out;
}

//convert inverse of homogeneous transform to Plucker transform
//HT=[R|t]
//P=
//[ R'          0 ]
//[-R'*skew3(t) R']
inline Mat6 invHTToPlucker(const HTransform& HT)
{
  Mat6 out;
  Mat3 Rt = HT.R.transpose();
  out.topLeftCorner<3,3>() = Rt;
  out.bottomLeftCorner<3,3>().noalias() = -Rt*skew(HT.t);
  out.topRightCorner<3,3>().setZero();
  out.bottomRightCorner<3,3>() = Rt;
  return out;
}

//FOR DEBUGGING
inline HTransform PluckerToHT(const Mat6& P)
{
  HTransform out;
  out.R = P.topLeftCorner<3,3>();
  out.t = unskew(P.bottomLeftCorner<3,3>()*out.R.transpose());
  return out;
}

//functions for multiplication with Plucker transforms
//can speed up computation given that Plucker B2 is always zeros

//multiply spatial vector by Plucker transform
//w=P*v
//[w0] = [P0 0 ] [v0]
//[w1]   [P1 P3] [v1]
//w0 = P0*v0
//w1 = P1*v0 + P3*v1
inline Vec6 multPluckerVec(const Mat6& P, const Vec6& v)
{
  Vec6 out;
  out.topRows<3>().noalias() = P.topLeftCorner<3,3>()*v.topRows<3>();
  out.bottomRows<3>().noalias() = P.bottomLeftCorner<3,3>()*v.topRows<3>() +
      P.bottomRightCorner<3,3>()*v.bottomRows<3>();
  return out;
}

//multiply spatial vector by transpose of Plucker transform
//w=P'*v
//[w0] = [P0' P1'] [v0]
//[w1]   [0   P3'] [v1]
//w0 = P0'*v0 + P1'*v1
//w1 = P3'*v1
inline Vec6 multPluckerTVec(const Mat6& P, const Vec6& v)
{
  Vec6 out;
  out.topRows<3>().noalias() = P.topLeftCorner<3,3>().transpose()*v.topRows<3>() +
      P.bottomLeftCorner<3,3>().transpose()*v.bottomRows<3>();
  out.bottomRows<3>().noalias() = P.bottomRightCorner<3,3>().transpose()*v.bottomRows<3>();
  return out;
}

//m is the mass
//c is the position of center of mass relative to reference frame
//I is the moment of inertia about the center of mass
//Is = [I + m*(C*C'), m*C ]
//     [m*C',         m*Id]
//where C = skew(c) and Id is identity matrix
inline Mat6 toSpatialInertia(const Real m, const Vec3& c, const Mat3& I)
{
  Mat3 C = skew(c);
  Mat6 out;
  out.topLeftCorner<3,3>().noalias() = I + m*(C*C.transpose());
  out.bottomLeftCorner<3,3>() = m*C;
  out.topRightCorner<3,3>() = m*C.transpose();
  out.bottomRightCorner<3,3>() = m*Mat3::Identity();
  return out;
}

void fromSpatialInertia(const Mat6& Is, Real& m, Vec3& c, Mat3& I)
{
  m = Is(3,3);
  Mat3 C = Is.topRightCorner<3,3>()/m;
  c = unskew(C);
  I = Is.topLeftCorner<3,3>() - m*(C*C.transpose());
}

//transform the coordinates of spatial inertia. returns:
//R = P'*I*P
//Plucker transforms & spatial inertia have special structure
//code generated using MATLAB symbolic toolbox
Mat6 multPluckerTInertiaPlucker(const Mat6& P, const Mat6& I);

} //namespace

#endif //_WMRDE_SPATIAL_H_


