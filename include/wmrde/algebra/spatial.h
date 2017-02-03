//utility for spatial vector algebra for rigid body dynamics
//Reference: Roy Featherstone, Rigid Body Dynamics Algorithms
//http://royfeatherstone.org/spatial/index.html#spatial-software

#ifndef _WMRDE_SPATIAL_H_
#define _WMRDE_SPATIAL_H_

#include <wmrde/algebra/transform.h>

namespace wmrde
{

class Vec6b //block representation of 6x1 vector
{
public:
  Vec3 b0, b1;

  Vec6b() {}
  Vec6b(
      const Vec3& block0,
      const Vec3& block1) :
        b0(block0),
        b1(block1) {}

  Vec6b operator+(const Vec6b& other) const
  {
    return Vec6b(
        b0 + other.b0,
        b1 + other.b1);
  }

  //FOR DEBUGGING
  inline Eigen::Matrix<Real,6,1> to6x1() const
  {
    Eigen::Matrix<Real,6,1> out;
    out.block(0,0,3,1) = b0;
    out.block(3,0,3,1) = b1;
    return out;
  }

  friend std::ostream &operator<<( std::ostream &output, const Vec6b &v)
  {
    output << v.to6x1();
    return output;
  }

  inline bool isApprox(const Vec6b& other, const Real prec = 0) const
  {
    return this->to6x1().isApprox(other.to6x1(), prec);
  }
};

class Mat6b //block representation of 6x6 matrix
{
public:
  //blocks in column major order:
  //M = [B0 B2]
  //    [B1 B3]

  Mat3 B0, B1, B2, B3;

  Mat6b() {}
  Mat6b(
      const Mat3& block0,
      const Mat3& block1,
      const Mat3& block2,
      const Mat3& block3) :
        B0(block0),
        B1(block1),
        B2(block2),
        B3(block3) {}

  Mat6b operator+(const Mat6b& other) const
  {
    return Mat6b(
        B0 + other.B0,
        B1 + other.B1,
        B2 + other.B2,
        B3 + other.B3);
  }

  inline Vec6b getColumn(const int idx) const
  {
    Vec6b out;
    if (idx < 3)
    {
      out.b0 << B0(0,idx), B0(1,idx), B0(2,idx);
      out.b1 << B1(0,idx), B1(1,idx), B1(2,idx);
    }
    else if (idx >= 3)
    {
      int idx_ = idx-3;
      out.b0 << B2(0,idx_), B2(1,idx_), B2(2,idx_);
      out.b1 << B3(0,idx_), B3(1,idx_), B3(2,idx_);
    }
    else
    {
      //TODO, handle out of bounds
    }
    return out;
  }

  //FOR DEBUGGING
  inline Eigen::Matrix<Real,6,6> to6x6() const
  {
    Eigen::Matrix<Real,6,6> out;
    out.block(0,0,3,3) = B0;
    out.block(3,0,3,3) = B1;
    out.block(0,3,3,3) = B2;
    out.block(3,3,3,3) = B3;
    return out;
  }

  friend std::ostream &operator<<( std::ostream &output, const Mat6b &M)
  {
    output << M.to6x6();
    return output;
  }

  inline bool isApprox(const Mat6b& other, const Real prec = 0) const
  {
    return this->to6x6().isApprox(other.to6x6(), prec);
  }
};

//CROSS PRODUCT OPERAITONS
//cross product with motion vector (m)
//d = v x m
//[d0] = [skew(v0) 0       ] [m0]
//[d1]   [skew(v1) skew(v0)] [m1]
//d0 = v0 x m0
//d1 = v1 x m0 + v0 x m1
inline Vec6b crossVec6bMotion(
    const Vec6b& v, const Vec6b& m)
{
  return Vec6b(
      v.b0.cross(m.b0),
      v.b1.cross(m.b0) + v.b0.cross(m.b1));
}

//cross product with force vector (f)
//d = v x f
//[d0] = [skew(v0) skew(v1)] [f0]
//[d1]   [0        skew(v0)] [f1]
//d0 = v0 x f0 + v1 x f1
//d1 = v0 x f1
inline Vec6b crossVec6bForce(
    const Vec6b& v, const Vec6b& f)
{
  return Vec6b(
      v.b0.cross(f.b0) + v.b1.cross(f.b1),
      v.b0.cross(f.b1));
}

//convert homogeneous transform to Plucker transform
//HT=[R|t]
//P=
//[R         0]
//[skew(t)*R R]
inline Mat6b HTToPlucker(const HTransform& HT)
{
  return Mat6b(HT.R, skew(HT.t)*HT.R, Mat3::Zero(), HT.R);
}

//convert inverse of homogeneous transform to Plucker transform
//HT=[R|t]
//P=
//[ R'          0 ]
//[-R'*skew3(t) R']
inline Mat6b invHTToPlucker(const HTransform& HT)
{
  Mat3 Rt = HT.R.transpose();
  return Mat6b(Rt, -Rt*skew(HT.t), Mat3::Zero(), Rt);
}

//FOR DEBUGGING
inline HTransform PluckerToHT(const Mat6b& P)
{
  return HTransform(P.B0, unskew(P.B1*P.B0.transpose()));
}

//functions for multiplication with Plucker transforms
//can speed up computation given that Plucker B2 is always zeros

//multiply spatial vector by Plucker transform
//w=P*v
//[w0] = [P0 0 ] [v0]
//[w1]   [P1 P3] [v1]
//w0 = P0*v0
//w1 = P1*v0 + P3*v1
inline Vec6b multPluckerVec6b(const Mat6b& P, const Vec6b& v)
{
  return Vec6b(P.B0*v.b0, P.B1*v.b0 + P.B3*v.b1);
}

//multiply spatial vector by transpose of Plucker transform
//w=P'*v
//[w0] = [P0' P1'] [v0]
//[w1]   [0   P3'] [v1]
//w0 = P0'*v0 + P1'*v1
//w1 = P3'*v1
inline Vec6b multPluckerTVec6b(const Mat6b& P, const Vec6b& v)
{
  return Vec6b(
      P.B0.transpose()*v.b0 + P.B1.transpose()*v.b1,
      P.B3.transpose()*v.b1);
}


//multiply spatial vector by Mat6b
//w = M*v
//[w0] = [M0 M2] [v0]
//[w1]   [M1 M3] [v1]
//w0 = M0*v0 + M2*v1
//w1 = M1*v0 + M3*v1
inline Vec6b multMatVec6b(const Mat6b& M, const Vec6b& v)
{
  return Vec6b(
      M.B0*v.b0 + M.B2*v.b1,
      M.B1*v.b0 + M.B3*v.b1);
}


//multiply matrix by transpose of Plucker transform
//R=P'*M
//[R0 R2] = [P0' P1'] [M0 M2]
//[R1 R3]   [0   P3'] [M1 M3]
//R0 = P0'*M0 + P1'*M1
//R1 = P3'*M1;
//R2 = P0'*M2 + P1'*M3
//R3 = P3'*M3
inline Mat6b multPluckerTMat6b(const Mat6b& P, const Mat6b& M)
{
  return Mat6b(
      P.B0.transpose()*M.B0 + P.B1.transpose()*M.B1,
      P.B3.transpose()*M.B1,
      P.B0.transpose()*M.B2 + P.B1.transpose()*M.B3,
      P.B3.transpose()*M.B3);

  //TODO, is this faster?
//  Mat6b out;
//  out.B0.noalias() = P.B0.transpose()*M.B0 + P.B1.transpose()*M.B1;
//  out.B1.noalias() = P.B3.transpose()*M.B1;
//  out.B2.noalias() = P.B0.transpose()*M.B2 + P.B1.transpose()*M.B3;
//  out.B3.noalias() = P.B3.transpose()*M.B3;
//  return out;
}


//multiply Plucker transform by matrix
//R=M*P
//[R0 R2] = [M0 M2] [P0 0 ]
//[R1 R3]   [M1 M3] [P1 P3]
//R0 = M0*P0 + M2*P1;
//R1 = M1*P0 + M3*P1;
//R2 = M2*P3;
//R3 = M3*P3;
inline Mat6b multMat6bPlucker(const Mat6b& M, const Mat6b& P)
{
  return Mat6b(
      M.B0*P.B0 + M.B2*P.B1,
      M.B1*P.B0 + M.B3*P.B1,
      M.B2*P.B3,
      M.B3*P.B3);

  //TODO, is this faster?
//  Mat6b out;
//  out.B0.noalias() = M.B0*P.B0 + M.B2*P.B1;
//  out.B1.noalias() = M.B1*P.B0 + M.B3*P.B1;
//  out.B2.noalias() = M.B2*P.B3;
//  out.B3.noalias() = M.B3*P.B3;
//  return out;
}

//TODO, is this faster?
inline void multMat6bPlucker(const Mat6b& M, const Mat6b& P, Mat6b& R)
{
  R.B0.noalias() = M.B0*P.B0 + M.B2*P.B1;
  R.B1.noalias() = M.B1*P.B0 + M.B3*P.B1;
  R.B2.noalias() = M.B2*P.B3;
  R.B3.noalias() = M.B3*P.B3;
}

//m is the mass
//c is the position of center of mass relative to reference frame
//I is the moment of inertia about the center of mass
//Is = [I + m*(C*C'), m*C ]
//     [m*C',         m*Id]
//where C = skew(c) and Id is identity matrix
inline Mat6b toSpatialInertia(const Real m, const Vec3& c, const Mat3& I)
{
  Mat3 C = skew(c);
  return Mat6b(I + m*(C*C.transpose()), m*C, m*C.transpose(), m*Mat3::Identity());
}

void fromSpatialInertia(const Mat6b& Is, Real& m, Vec3& c, Mat3& I)
{
  m = Is.B3(0,0);
  Mat3 C = Is.B2/m;
  c = unskew(C);
  I = Is.B0 - m*(C*C.transpose());
}


} //namespace

#endif //_WMRDE_SPATIAL_H_


