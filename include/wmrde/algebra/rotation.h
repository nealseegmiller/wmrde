//utility for converting between rotation representations

#ifndef _WMRDE_ROTATION_H_
#define _WMRDE_ROTATION_H_

#include <wmrde/algebra/linalg3.h>
#include <Eigen/Geometry>

namespace wmrde
{
typedef Eigen::Quaternion<Real> Quaternion;
typedef Eigen::AngleAxis<Real> AngleAxis;

inline Real radToDeg(Real rad) { return rad*180.0/M_PI; }
inline Real degToRad(Real deg) { return deg*M_PI/180.0; }

/*!
 * compute 3D rotation matrix about x-axis
 */
inline Mat3 Rotx(const Real angle)
{
  Real ca = cos(angle);
  Real sa = sin(angle);
  Mat3 out;
  out << 1, 0,   0,
         0, ca,-sa,
         0, sa, ca;
  return out;
}

/*!
 * compute 3D rotation matrix about y-axis
 */
inline Mat3 Roty(const Real angle)
{
  Real ca = cos(angle);
  Real sa = sin(angle);
  Mat3 out;
  out << ca, 0, sa,
         0,  1, 0,
        -sa, 0, ca;
  return out;
}

/*!
 * compute 3D rotation matrix about z-axis
 */
inline Mat3 Rotz(const Real angle)
{
  Real ca = cos(angle);
  Real sa = sin(angle);
  Mat3 out;
  out << ca, -sa, 0,
         sa,  ca, 0,
         0,   0,  1;
  return out;
}

/*!
 * convert Euler angles (in radians) to rotation matrix. return value is equivalent to:
 * R = Rotz(yaw)*Roty(pitch)*Rotx(roll)
 * but is computed more efficiently
 */
inline Mat3 eulerToRot(const Real roll, const Real pitch, const Real yaw)
{
  Real sR=sin(roll);  Real cR=cos(roll);
  Real sP=sin(pitch); Real cP=cos(pitch);
  Real sY=sin(yaw);   Real cY=cos(yaw);

  Mat3 out;

  out(0,0) = cP*cY;
  out(0,1) = cY*sP*sR - cR*sY;
  out(0,2) = sR*sY + cR*cY*sP;
  out(1,0) = cP*sY;
  out(1,1) = cR*cY + sP*sR*sY;
  out(1,2) = cR*sP*sY - cY*sR;
  out(2,0) =-sP;
  out(2,1) = cP*sR;
  out(2,2) = cP*cR;

  //TODO, why is this slow?
//  out << cP*cY, cY*sP*sR - cR*sY, sR*sY + cR*cY*sP,
//         cP*sY, cR*cY + sP*sR*sY, cR*sP*sY - cY*sR,
//        -sP,    cP*sR,            cP*cR;

  return out;
}

/*!
 * convert rotation matrix to Euler angles
 */
inline void rotToEuler(const Mat3 R,
    Real& roll, Real& pitch, Real& yaw) //output
{
  yaw = atan2(R(1,0), R(0,0));
  roll = atan2(R(2,1), R(2,2));
  Real sR=sin(roll);
  Real cR=cos(roll);
  Real cP;
  if (fabs(cR) > fabs(sR)) {
    cP = R(2,2)/cR; //singularity at rol = +/- pi/2
  } else {
    cP = R(2,1)/sR; //singularity at rol = 0,+/-pi
  }
  pitch = atan2(-R(2,0), cP); //pit
}


/*
 * An alternative Euler angle to rotation implementation using AngleAxis
 * is suggested in the Eigen documentation:
 * https://eigen.tuxfamily.org/dox/classEigen_1_1AngleAxis.html
 * use this implementation for unit testing and benchmarking
 */
inline Mat3 RotxTest(const Real angle) { return AngleAxis(angle, Vec3::UnitX()).toRotationMatrix(); }
inline Mat3 RotyTest(const Real angle) { return AngleAxis(angle, Vec3::UnitY()).toRotationMatrix(); }
inline Mat3 RotzTest(const Real angle) { return AngleAxis(angle, Vec3::UnitZ()).toRotationMatrix(); }
inline Mat3 eulerToRotTest(const Real roll, const Real pitch, const Real yaw)
{
  return (
      AngleAxis(yaw, Vec3::UnitZ())*
      AngleAxis(pitch, Vec3::UnitY())*
      AngleAxis(roll, Vec3::UnitX()) ).toRotationMatrix();
}

/*!
 * transform from angular velocity (in body coords) to Euler angle rates
 * \param roll roll angle
 * \param pitch pitch angle
 * \return transform T such that (euler angle rates) = T*(angular velocity)
 */
inline Mat3 velToEulerrateTransform(const Real roll, const Real pitch)
{
  Real sR = sin(roll);
  Real cR = cos(roll);
  Real cP = cos(pitch);
  Real tP = tan(pitch);

  Mat3 T;
  T << 1, sR*tP, cR*tP,
       0, cR,   -sR,
       0, sR/cP, cR/cP;
  return T;
}

/*!
 * transform from Euler angle rates to angular velocity (in body coords)
 * \param roll roll angle
 * \param pitch pitch angle
 * \return transform T such that (angular velocity) = T*(euler angle rates)
 */
inline Mat3 eulerrateToVelTransform(const Real roll, const Real pitch)
{
  Real sR = sin(roll);
  Real cR = cos(roll);
  Real sP = sin(pitch);
  Real cP = cos(pitch);

  Mat3 T;
  T << 1, 0, -sP,
       0, cR, sR*cP,
       0,-sR, cR*cP;
  return T;
}

/*!
 * Step the rotation of an object wrt world (represented as quaternion)
 * given angular velocity and time step.
 * \param q rotation quaternion before step
 * \param angular_vel angular velocity *in world coords*
 * \param dt time step size
 * \return rotation quaternion after step
 */
inline Quaternion stepRotationQuat(
    const Quaternion q,
    const Vec3& angular_vel,
    const Real dt)
{
  Real norm = angular_vel.norm();
  if (norm > 0)
  {
    Quaternion dq(AngleAxis(norm*dt, angular_vel/norm));
    return (dq*q).normalized();
  }
  return q;
}

} //namespace

#endif //_WMRDE_TRANSFORM_H_
