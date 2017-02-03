#include <gtest/gtest.h>

#include <wmrde/algebra/rotation.h>
#include <wmrde/algebra/spatial.h>
#include <wmrde/util/string_format.h>

//uncomment one to enable/disable print output
#define DEBUG_INFO(x) do { std::cout << x << std::endl; } while(0)
//#define DEBUG_INFO(x) do { } while(0)

using namespace wmrde;

Real epsilon = std::numeric_limits<Real>::epsilon()*2;
inline bool approx(const Real a, const Real b) { return std::abs(a-b) < epsilon; }

TEST(TestSuite, rotation)
{
  Mat3 Rx,Ry,Rz,R;
  Real roll = degToRad(10.0);
  Real pitch = degToRad(20.0);
  Real yaw = degToRad(30.0);

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
  HTransform T( eulerToRot(0.0,0.0,degToRad(10.0)), Vec3(1.0,2.0,3.0) );
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

//TODO, why does including the spatial test suite cause
// some benchmark times to increase in previous test suites?!
// affects benchmarks even if DISABLED_. need to comment it out!
#define RUN_SPATIAL_TEST 1
#if RUN_SPATIAL_TEST

typedef Eigen::Matrix<Real,6,6> Mat6;
typedef Eigen::Matrix<Real,6,1> Vec6;

TEST(TestSuite, spatial) //prepend DISABLED_ to name to disable
{
  //validate functions to construct Plucker transform
  //from homogeneous transform
  HTransform HT( eulerToRot(0,0,degToRad(10.0)), Vec3(1.0,2.0,3.0) );
  Mat6b P = HTToPlucker(HT);
  Mat6b Pinv = invHTToPlucker(HT);
  HTransform HT_ = PluckerToHT(P);

  DEBUG_INFO("HT = \n" << HT);
  DEBUG_INFO("P(HT) = \n" << P);
  DEBUG_INFO("P(inverse HT) = \n" << Pinv);
  DEBUG_INFO("HT(P) = \n" << HT_);

  Mat6 P_nb = P.to6x6(); //non-block
  Mat6 Id = Pinv.to6x6()*P_nb; //expect identity

  DEBUG_INFO("inverse(P)*P = \n" << Id);

  EXPECT_TRUE( Id.isApprox(Mat6::Identity()) );
  EXPECT_TRUE( HT_.isApprox(HT) );

  //validate getColumn
  for (size_t i = 0; i < 6; i++)
  {
    Vec6b col = P.getColumn(i);
    EXPECT_TRUE(P_nb.col(i).isApprox(col.to6x1()));
//    std::cout << "P(:," << i << ") = \n" << col << std::endl;
  }

  {
    //validate Plucker-vector multiplication
    Vec6b v( Vec3(1.0,2.0,3.0), Vec3(4.0,5.0,6.0) );
    Vec6b Pv = multPluckerVec6b(P,v);
    Vec6b PTv = multPluckerTVec6b(P,v);

    DEBUG_INFO("v = \n" << v);
    DEBUG_INFO("P*v = \n" << Pv);
    DEBUG_INFO("P'*v = \n" << PTv);

    Vec6 v_nb = v.to6x1(); //non-block
    EXPECT_TRUE(Pv.to6x1().isApprox( P_nb*v_nb ));
    EXPECT_TRUE(PTv.to6x1().isApprox( P_nb.transpose()*v_nb ));
  }

  {
    //validate Plucker-matrix multiplication
    Real mass = 10.0;
    Mat6b I = toSpatialInertia(mass, Vec3(1.0, 2.0, 3.0), Mat3::Identity());

    Mat6b PTI = multPluckerTMat6b(P,I);
    Mat6b IP = multMat6bPlucker(I,P);

    DEBUG_INFO("I = \n" << I);
    DEBUG_INFO("P'*I = \n" << PTI);
    DEBUG_INFO("I*P = \n" << IP);

    Mat6 I_nb = I.to6x6();
    EXPECT_TRUE(PTI.to6x6().isApprox(P_nb.transpose()*I_nb));
    EXPECT_TRUE(IP.to6x6().isApprox(I_nb*P_nb));
  }
}
#endif //RUN_SPATIAL_TEST

/*
TEST(TestSuite, matrix)
{
  printf("benchmark matrix.h");

  const int SIZE = 16;

  //matrix.h matrices
  Real A[SIZE*SIZE];
  Real B[SIZE*SIZE];
  Real B0[SIZE*SIZE]; //backup, to initialize Eigen
  Real B2[SIZE*SIZE];
  Real b[SIZE];
  Real r[SIZE];

  for (int i=0; i<SIZE*SIZE; i++)
  {
    A[i] = (Real) i+1;
    B[i] = 1;
    B0[i] = B[i];
  }

  for (int i=0; i<SIZE; i++) { b[i] = 1; }

//  std::cout << "A=\n"; printMatReal(SIZE,SIZE,A,-1,-1);
//  std::cout << "B=\n"; printMatReal(SIZE,SIZE,B,-1,-1);
//  std::cout << "b=\n"; printMatReal(SIZE,1,b,-1,-1);

  //for block operations
  int ri = 2;
  int ci = 3;
  int brows = SIZE/2;
  int bcols = SIZE/2;

  Real val = 2.0;
  Real eps = 1e-6;
//  Real m = 1.0 + eps;
  Real m = 1.0;


  size_t num_iter = 1e6;
  Timer timer;

  {
    timer.start();
    for (size_t i=0; i<num_iter; i++)
    {
//      setMat(SIZE,SIZE,val,B);
      //setMatRow(SIZE,SIZE,ri,val,B);
      //setMatCol(SIZE,ci,val,B);
      //setMatBlock(SIZE,ri,ci,brows,bcols,val,B);

  //    mulcMat(SIZE,SIZE,m,B);
      //mulcMatRow(SIZE,SIZE,ri,m,B);
      //mulcMatCol(SIZE,ci,m,B);
      //mulcMatBlock(SIZE,ri,ci,brows,bcols,m,B);

      //copyMat(SIZE,SIZE,A,B);
      //copyMatRow(SIZE,SIZE,ri,A,SIZE,ri,B);
      //copyMatCol(SIZE,ci,A,ci,B);
      //copyMatBlock(SIZE,ri,ci,brows,bcols,A, SIZE,ri,ci,B);

      //copyTMat(SIZE,SIZE,A,B);

//      addmMat(SIZE,SIZE,A,m,B);
      //addmMatRow(SIZE,SIZE,ri,A,SIZE,ri,m,B);
      //addmMatCol(SIZE,ci,A,ci,m,B);
  //    addmMatBlock(SIZE,ri,ci,brows,bcols,A, SIZE,ri,ci,m,B);

//      multMatVec(SIZE,SIZE,A,b,m,r);
      multMatTVec(SIZE,SIZE,A,b,m,B);
//      multMatMat(SIZE,SIZE,A,SIZE,A,m,B);
  //    multMatTMat(SIZE,SIZE,A,SIZE,A,m,B);
    }
    timer.stop();

//    std::cout << "(matrix.h) B=\n"; printMatReal(SIZE,SIZE,B,0,-1);
    std::cout << "(matrix.h) r=\n"; printMatReal(SIZE,1,r,0,-1);
    printf("for %zu iterations using matrix.h, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

  {
    //Eigen

    //dynamic Eigen matrices
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> A_(SIZE,SIZE);
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> B_(SIZE,SIZE);
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> b_(SIZE,1);
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> r_(SIZE,1);

    //fixed Eigen matrices
//    Eigen::Matrix<Real,SIZE,SIZE> A_;
//    Eigen::Matrix<Real,SIZE,SIZE> B_;
//    Eigen::Matrix<Real,SIZE,1> b_;
//    Eigen::Matrix<Real,SIZE,1> r_;

    memcpy(A_.data(), A, sizeof(Real)*SIZE*SIZE);
    memcpy(B_.data(), B0, sizeof(Real)*SIZE*SIZE);
    memcpy(b_.data(), b, sizeof(Real)*SIZE);

    //std::cout << "Eigen:\n";
    //std::cout << "A_=\n" << A_ << std::endl;
    //std::cout << "B_=\n" << B_ << std::endl;
    //std::cout << "b_=\n" << b_ << std::endl;

    timer.start();
    for (size_t i=0; i<num_iter; i++)
    {
//      B_.setConstant(val);
      //B_.row(ri).setConstant(val);
      //B_.col(ci).setConstant(val);
      //B_.block(ri,ci,brows,bcols).setConstant(val);

//      B_ *= m;
      //B_.row(ri) *= m;
      //B_.col(ci) *= m;
      //B_.block(ri,ci,brows,bcols) *= m;

      //B_ = A_;
      //B_.row(ri) = A_.row(ri);
      //B_.col(ci) = A_.col(ci);
//      B_.block(ri,ci,brows,bcols) = A_.block(ri,ci,brows,bcols);

      //B_ = A_.transpose();

//      B_ += (A_ * m);
      //B_.row(ri) += (A_.row(ri) * m);
      //B_.col(ci) += (A_.col(ci) * m);
      //B_.block(ri,ci,brows,bcols) += (A_.block(ri,ci,brows,bcols) * m);

//      r_.noalias() = A_*b_;
      r_.noalias() = A_.transpose()*b_;
//      B_.noalias = A_*A_;
//      B_.noalias() = A_.transpose()*A_;
    }
    timer.stop();

//    std::cout << "(Eigen) B=\n" << B_ << std::endl;
    std::cout << "(Eigen) r=\n" << r_ << std::endl;
    printf("for %zu iterations using Eigen, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }
}
*/

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
