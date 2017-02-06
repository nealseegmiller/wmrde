#include <gtest/gtest.h>

#include <wmrde/algebra/rotation.h>
#include <wmrde/algebra/spatial.h>
#include <wmrde/algebra/spatial_block.h>
#include <wmrde/algebra/random.h>
#include <wmrde/util/string_format.h>

//uncomment one to enable/disable print output
#define DEBUG_INFO(x) do { std::cout << x << std::endl; } while(0)
//#define DEBUG_INFO(x) do { } while(0)

using namespace wmrde;

Real epsilon = std::numeric_limits<Real>::epsilon()*10.0;
inline bool approx(const Real a, const Real b)
{
  return std::abs(a-b) < epsilon;
}


TEST(TestSuite, rotation)
{
  Mat3 Rx,Ry,Rz,R;
  Real roll = randAngle();
  Real pitch = randAngle();
  Real yaw = randAngle();

  //calculate rotation matrices from Euler angles
  Rx = Rotx(roll);
  Ry = Roty(pitch);
  Rz = Rotz(yaw);
  R = eulerToRot(roll,pitch,yaw);

  DEBUG_INFO(string_format("Rotx(%f) = \n", roll) << Rx);
  DEBUG_INFO(string_format("Roty(%f) = \n", pitch) << Ry);
  DEBUG_INFO(string_format("Rotz(%f) = \n", yaw) << Rz);
  DEBUG_INFO(string_format("eulerToRot(%f,%f,%f) = \n", roll, pitch, yaw) << R);

  //to validate, calculate rotation matrices from Euler angles using AngleAxis
  //as suggested in the Eigen documentation:
  //https://eigen.tuxfamily.org/dox/classEigen_1_1AngleAxis.html
  Mat3 Rx2 = RotxTest(roll);
  Mat3 Ry2 = RotyTest(pitch);
  Mat3 Rz2 = RotzTest(yaw);
  Mat3 R2 = eulerToRotTest(roll,pitch,yaw);

  EXPECT_TRUE( Rx.isApprox(Rx2) );
  EXPECT_TRUE( Ry.isApprox(Ry2) );
  EXPECT_TRUE( Rz.isApprox(Rz2) );
  EXPECT_TRUE( R.isApprox(R2) );

  //validate rotToEuler()
  Real roll_, pitch_, yaw_;
  rotToEuler(R2, roll_, pitch_, yaw_);
  DEBUG_INFO(string_format("rotToEuler(): %f, %f, %f", roll_, pitch_, yaw_));

  EXPECT_TRUE(
      approx(roll,roll_) &&
      approx(pitch,pitch_) &&
      approx(yaw,yaw_));
}

TEST(TestSuite, transform)
{
  //validate composeWith()
  HTransform T = randomHTransform();
  HTransform TT = T.composeWith(T);

  DEBUG_INFO("T = \n" << T);
  DEBUG_INFO("T*T = \n" << TT);

  //validate using 4x4 matrix multiplication
  Eigen::Matrix<Real,4,4> T_nb = T.to4x4(); //non-block
  EXPECT_TRUE( TT.to4x4().isApprox(T_nb*T_nb) );

  //validate composeInvWith()
  HTransform Tid = T.composeInvWith(T);
  HTransform Tid2; Tid2.setIdentity();

  DEBUG_INFO("inverse(T)*T = \n" << Tid);

  EXPECT_TRUE( Tid.isApprox(Tid2, epsilon) );
  EXPECT_TRUE( Tid.isApprox(T.inverse().composeWith(T)) );

  //validate applyTo() and applyInvTo()
  EXPECT_TRUE( T.applyTo(T.t).isApprox(TT.t) );
  EXPECT_TRUE( T.applyInvTo(T.t).isApprox(Vec3::Zero()) );
}


TEST(TestSuite, spatial_convert)
{
  //validate functions to construct Plucker transform
  //from homogeneous transform
  HTransform HT = randomHTransform();
  Mat6 P = HTToPlucker(HT);
  Mat6 Pinv = invHTToPlucker(HT);
  HTransform HT_ = PluckerToHT(P);

  DEBUG_INFO("HT = \n" << HT);
  DEBUG_INFO("P(HT) = \n" << P);
  DEBUG_INFO("P(inverse HT) = \n" << Pinv);
  DEBUG_INFO("HT(P) = \n" << HT_);

  Mat6 Pinv_P = Pinv*P;
  DEBUG_INFO("inverse(P)*P = \n" << Pinv_P);

  EXPECT_TRUE( Pinv_P.isApprox(Mat6::Identity()) );
  EXPECT_TRUE( HT_.isApprox(HT) );
}

TEST(TestSuite, spatial_cross)
{
  //validate cross products
  Vec6 v; v << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  Vec6 m; m << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0;
  Vec6 f = m;

  Vec6 vxm = crossVecMotion(v,m);
  Vec6 vxf = crossVecForce(v,f);

  DEBUG_INFO("v = \n" << v);
  DEBUG_INFO("m = \n" << m);
  DEBUG_INFO("f = \n" << f);
  DEBUG_INFO("v x m = \n" << vxm);
  DEBUG_INFO("v x f = \n" << vxf);

  Vec6 vxm_expect; vxm_expect << -1.0, 2.0, -1.0, -2.0, 4.0, -2.0; //calculated by hand
  Vec6 vxf_expect; vxf_expect << -2.0, 4.0, -2.0, -4.0, 8.0, -4.0; //calculated by hand

  EXPECT_TRUE( vxm.isApprox(vxm_expect) );
  EXPECT_TRUE( vxf.isApprox(vxf_expect) );
}

TEST(TestSuite, spatial_plucker_vector)
{
  //validate Plucker-vector multiplication
  Mat6 P = randomPlucker();
  Vec6 v; v << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  Vec6 Pv = multPluckerVec(P,v);
  Vec6 PTv = multPluckerTVec(P,v);

  DEBUG_INFO("P = \n" << P);
  DEBUG_INFO("v = \n" << v);
  DEBUG_INFO("P*v = \n" << Pv);
  DEBUG_INFO("P'*v = \n" << PTv);

  EXPECT_TRUE(Pv.isApprox(P*v));
  EXPECT_TRUE(PTv.isApprox(P.transpose()*v));
}

//validate P'*I*P calculation using multPluckerTInertiaPlucker(),
//which uses optimized code generated with symbolic math toolbox.
//TODO, eliminate multPluckerTInertiaPlucker function and this test because
//benchmarking shows it's not much faster
TEST(TestSuite, spatial_plucker_inertia)
{
  //validate Plucker-matrix multiplication
  Mat6 P = randomPlucker();
  Mat6 I = randomInertia();

  Mat6 PtIP = multPluckerTInertiaPlucker(P,I);
  Mat6 PtIP_expect = P.transpose()*I*P;

  DEBUG_INFO("P = \n" << P);
  DEBUG_INFO("I = \n" << I);
  DEBUG_INFO("P'*I*P = \n" << PtIP);

  EXPECT_TRUE( PtIP.isApprox(PtIP_expect,epsilon) );
}

//validate P'*I*P calculation using deprecated Mat6b class
//in spatial_block.h
//TODO, eliminate spatial_block.h and this test because
//benchmarking shows it's not much faster.
TEST(TestSuite, spatial_plucker_inertia_block)
{
  //validate Plucker-matrix multiplication
  Mat6b P; P.from6x6(randomPlucker());
  Mat6b I; I.from6x6(randomInertia());

  Mat6b M = multMatPlucker6b(I,P);
  Mat6b PtIP = multPluckerTMat6b(P,M);

  Mat6 P_nb = P.to6x6();
  Mat6 I_nb = I.to6x6();
  Mat6 PtIP_expect = P_nb.transpose()*I_nb*P_nb;

  DEBUG_INFO("P = \n" << P);
  DEBUG_INFO("I = \n" << I);
  DEBUG_INFO("P'*I*P = \n" << PtIP);

  EXPECT_TRUE( PtIP.to6x6().isApprox(PtIP_expect,epsilon) );
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
