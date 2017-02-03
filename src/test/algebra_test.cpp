#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <wmrde/algebra/rotation.h>
#include <wmrde/algebra/spatial.h>
#include <wmrde/util/string_format.h>
#include <wmrde/util/timer.h>

using namespace wmrde;

Real epsilon = std::numeric_limits<Real>::epsilon()*2;
inline bool approx(const Real a, const Real b)
{
  return std::abs(a-b) < epsilon;
}

TEST(TestSuite, linalg3)
{
  Mat3 A, B, C;
  for (int i = 0; i < A.size(); i++)
  {
    A.data()[i] = i;
    B.data()[i] = i+1;
  }

  std::cout << "A = \n" << A << std::endl;
  std::cout << "B = \n" << B << std::endl;
  std::cout << "A*B = \n" << A*B << std::endl;
  std::cout << "A'*B = \n" << A.transpose()*B << std::endl;

//  EXPECT_TRUE(???);

  size_t num_iter = 1e7;
  Timer timer;

  //benchmark A*B
  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++) { C.noalias() = A*B; }
  timer.stop();
  printf("%zu iterations of Mat3 A*B took %f ms\n", num_iter, timer.elapsedTimeMs());

  //benchmark A'*B
  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++) { C.noalias() = A.transpose()*B; }
  timer.stop();
  printf("%zu iterations of Mat3 A'*B took %f ms\n", num_iter, timer.elapsedTimeMs());
}

inline Mat3 RotxTest(const Real angle) { return AngleAxis(angle, Vec3::UnitX()).toRotationMatrix(); }
inline Mat3 RotyTest(const Real angle) { return AngleAxis(angle, Vec3::UnitY()).toRotationMatrix(); }
inline Mat3 RotzTest(const Real angle) { return AngleAxis(angle, Vec3::UnitZ()).toRotationMatrix(); }
inline Mat3 RotTest(const Real roll, const Real pitch, const Real yaw)
{
  return (
      AngleAxis(yaw, Vec3::UnitZ())*
      AngleAxis(pitch, Vec3::UnitY())*
      AngleAxis(roll, Vec3::UnitX()) ).toRotationMatrix();
}

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

  std::cout << string_format("Rotx(%f) = \n", roll) << Rx << std::endl;
  std::cout << string_format("Roty(%f) = \n", pitch) << Ry << std::endl;
  std::cout << string_format("Rotz(%f) = \n", yaw) << Rz << std::endl;
  std::cout << string_format("eulerToRot(%f,%f,%f) = \n", roll, pitch, yaw) << R << std::endl;

  //to validate, calculate rotation matrices from Euler angles using AngleAxis
  //as suggested in the Eigen documentation:
  //https://eigen.tuxfamily.org/dox/classEigen_1_1AngleAxis.html
  Mat3 Rx2 = RotxTest(roll);
  Mat3 Ry2 = RotyTest(pitch);
  Mat3 Rz2 = RotzTest(yaw);
  Mat3 R2 = RotTest(roll,pitch,yaw);

  EXPECT_TRUE( Rx.isApprox(Rx2) ) << "Rotx() failed, expected:\n" << Rx;
  EXPECT_TRUE( Ry.isApprox(Ry2) ) << "Roty() failed, expected:\n" << Ry;
  EXPECT_TRUE( Rz.isApprox(Rz2) ) << "Rotz() failed, expected:\n" << Rz;
  EXPECT_TRUE( R.isApprox(R2) ) << "eulerToRot() failed, expected:\n" << R2;

  //validate rotToEuler()
  Real roll_, pitch_, yaw_;
  rotToEuler(R2, roll_, pitch_, yaw_);
  printf("rotToEuler() roll = %f, pitch = %f, yaw = %f\n", roll_, pitch_, yaw_);
  EXPECT_TRUE(approx(roll,roll_) && approx(pitch,pitch_) && approx(yaw,yaw_)) << "rotToEuler() failed.\n";

  //benchmark computation time of rotation.h functions vs. using AngleAxis
  size_t num_iter = 1e7;
  Timer timer;

  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    R = eulerToRot(roll,pitch,yaw);
//    Rx = Rotx(roll);
  }
  timer.stop();
  printf("%zu iterations of euler to Rotation took %f ms\n", num_iter, timer.elapsedTimeMs());

  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    R2 = RotTest(roll,pitch,yaw);
//    Rx2 = RotxTest(roll);
  }
  timer.stop();
  printf("%zu iterations of euler to Rotation using AngleAxis took %f ms\n", num_iter, timer.elapsedTimeMs());
}

TEST(TestSuite, transform)
{
  //validate composeWith()
  HTransform T( eulerToRot(0.0,0.0,degToRad(10.0)), Vec3(1.0,2.0,3.0) );
  HTransform TT = T.composeWith(T);

  std::cout << "T = \n" << T << std::endl;
  std::cout << "T*T = \n" << TT << std::endl;

  //validate using 4x4 matrix multiplication
  Eigen::Matrix<Real,4,4> T4x4 = T.to4x4();
  EXPECT_TRUE( TT.to4x4().isApprox(T4x4*T4x4) );

  //validate composeInvWith()
  HTransform Tid = T.composeInvWith(T);
  HTransform Tid2; Tid2.setIdentity();
  std::cout << "inv(T)*T = \n" << Tid << std::endl;
  EXPECT_TRUE( Tid.isApprox(Tid2, epsilon) );
  EXPECT_TRUE( Tid.isApprox(T.inverse().composeWith(T)) );

  //validate applyTo() and applyInvTo()
  EXPECT_TRUE( T.applyTo(T.t).isApprox(TT.t) );
  EXPECT_TRUE( T.applyInvTo(T.t).isApprox(Vec3::Zero()) );

  //benchmark composeWith() agains 4x4 matrix multiplication
  size_t num_iter = 1e7;
  Timer timer;

  //TODO, write a benchmarking macro?
  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    TT = T.composeWith(T);
  }
  timer.stop();
  printf("%zu iterations of composeWith() took %f ms\n", num_iter, timer.elapsedTimeMs());

  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    Tid = T.composeInvWith(T);
  }
  timer.stop();
  printf("%zu iterations of composeInvWith() took %f ms\n", num_iter, timer.elapsedTimeMs());

  Eigen::Matrix<Real,4,4> Mat4x4;
  timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    Mat4x4.noalias() = T4x4 * T4x4;
  }
  timer.stop();
  printf("%zu iterations of 4x4 matrix multiplication took %f ms\n", num_iter, timer.elapsedTimeMs());
}

/*
TEST(TestSuite, addMat3)
{
  printf("benchmark addMat3\n");
  Mat3 A,B,R,R2;
  Real eps = 1e-6;
  setMat3Diagonal(eps, eps, eps, A);
  setMat3Diagonal(eps, eps, eps, B);

  Mat3 A0; copyMat3(A,A0); //backup

  size_t num_iter = 1e7;
  Timer timer;
  {
    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      addMat3(A,B,A);
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of addMat3, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
    copyMat3(A,R); //first result
  }

  //verify with Eigen library
  {
    Eigen::Matrix<Real,3,3> A_, B_;
    copyMat3ToArray(A0, A_.data());
    copyMat3ToArray(B, B_.data());

    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      A_ = A_ + B_;
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of Eigen 3x3 matrix addition, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copyArrayToMat3(A_.data(), R2);
  }

  printf("R = \n%s", Mat3ToString(R).c_str());
  printf("R2 = \n%s", Mat3ToString(R2).c_str());

  EXPECT_TRUE(Mat3Equal(R,R2)); //second result
}

TEST(TestSuite, multMatMat3)
{
  Mat3 A,B,R,R2;
  Real eps = 1e-6;
  setMat3Diagonal(1+eps, 1+eps, 1+eps, A);
  setMat3Diagonal(1+eps, 1+eps, 1+eps, B); B[1] = eps;

  Mat3 A0; copyMat3(A,A0);

  size_t num_iter = 1e7;
  Timer timer;
  {
    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
//      multMatMat3(A,B,R);
      multMatTMat3(A,B,R);
      copyMat3(R,A);
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of multMatMat3, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

  //verify with Eigen library
  {
    copyMat3(A0,A);

    Eigen::Matrix<Real,3,3> A_, B_;
    copyMat3ToArray(A, A_.data());
    copyMat3ToArray(B, B_.data());

    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
//      A_ = A_*B_;
      A_ = A_.transpose()*B_;
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of Eigen 3x3 matrix multiplication, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copyArrayToMat3(A_.data(), R2);
  }

  printf("R = \n%s", Mat3ToString(R).c_str());
  printf("R2 = \n%s", Mat3ToString(R2).c_str());

  EXPECT_TRUE(Mat3Equal(R,R2));
}

TEST(TestSuite, composeHT)
{
  printf("benchmark composeHT\n");
  VecEuler euler;
  Vec3 translation;
  HomogeneousTransform HT,HTR;

  setEuler(DEGTORAD(0),DEGTORAD(0),DEGTORAD(0.1),euler);
  Real eps = 1e-6;
  setVec3(eps,eps,eps,translation);
  poseToHT(euler,translation,HT);

  HomogeneousTransform HT0; copyHT(HT, HT0); //backup

  std::cout << "HT=\n"; printHT(HT,-1,-1);

  //test compose HT
  int num_iter = (int) 1e7;
  Timer timer;
  {
    //time it
    timer.start();
    for (int i=0; i<num_iter; i++)
    {
      composeHT(HT,HT0,HTR);
      copyHT(HTR,HT);
    }
    timer.stop();

    printf("HT*HT = \n");
    printHT(HTR,-1,-1);

    printf("for %zu iterations using composeHT, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

  //validate using Eigen::Matrix<Real,3,3>
  {
    Eigen::Matrix<Real,3,3> R0_,R_; //rotation
    Eigen::Matrix<Real,3,1> t0_,t_; //translation

    copyMat3ToArray(HT0,R_.data());
    copy3(HT0+COL3,t_.data());
    std::cout << "Eigen HT=\n" << R_ << "\n" << t_ << std::endl;

    R0_ = R_; //backup
    t0_ = t_;

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
      t_ = R_*t0_ + t_;
      R_ = R_*R0_;
    }
    timer.stop();

    std::cout << "Eigen HT*HT=\n" << R_ << "\n" << t_ << std::endl;
    printf("for %zu iterations using Eigen 3x3, 3x1 blocks, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

  //validate using Eigen::Matrix<Real,4,4>
  {
    Eigen::Matrix<Real,4,4> HT_, HT0_;
    HT_.setZero();

    for (int i=0; i<4; i++) { copy3(HT0+(i*VEC3_SIZE), HT_.data() + (i*4)); }
    HT_(3,3) = 1.0;
    std::cout << "Eigen HT=\n" << HT_ << std::endl;
    HT0_ = HT_; //backup

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
      HT_ = HT_*HT0_;
    }
    timer.stop();

    std::cout << "Eigen HT*HT=\n" << HT_ << std::endl;
    printf("for %zu iterations using Eigen 4x4 matrices, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

}

TEST(TestSuite, multPluckerVec6b)
{
  printf("benchmark multPluckerVec6b");

  //construct Homogeneous Transform
  VecEuler euler;
  Vec3 translation;
  HomogeneousTransform HT;

  setEuler(DEGTORAD(10),DEGTORAD(20),DEGTORAD(30),euler);
  setVec3(1,2,3,translation);
  poseToHT(euler,translation,HT);

  std::cout << "HT=\n";
  printHT(HT,-1,-1);

  //convert to Plucker transform
  Mat6b P;
  HTToPlucker(HT,P);
  std::cout << "P=\n"; printMat6b(P,-1,-1);

  Vec6b m,v;
  setVec6b(1,2,3,4,5,6,m);
  std::cout << "m=\n"; printVec6b(m,-1,-1);

  Vec6b m0; copyVec6b(m,m0);

  size_t num_iter = 1e7;
  Timer timer;
  //multiply vector by Plucker transform
  {
    //time it

    timer.start();
    for (size_t i=0; i<num_iter; i++)
    {
      multPluckerVec6b(P,m,v);
      copyVec6b(v,m);
//      multPluckerTVec6b(P,m,v);
//      copyVec6b(v,m);
    }
    timer.stop();

    std::cout << "P*m=\n"; printVec6b(v,-1,-1);
//    std::cout << "P'*m=\n"; printVec6b(v,-1,-1);

    printf("for %zu iterations using spatial.h, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

  {
    //using Eigen 3x3 blocks
    Eigen::Matrix<Real,3,3> B0_, B1_, B3_;
    Eigen::Matrix<Real,3,1> m0_, m1_, m0_backup_, m1_backup_;

    //copy to Eigen matrices
    copyMat3ToArray(P,B0_.data());
    copyMat3ToArray(P+BLOCK1, B1_.data());
    copyMat3ToArray(P+BLOCK3, B3_.data());
    copy3(m0,m0_.data());
    copy3(m0+VEC3_SIZE, m1_.data());

    timer.start();
    for (size_t i=0; i<num_iter; i++)
    {
      m0_backup_ = m0_;
      m0_ = B0_*m0_;
      m1_ = B1_*m0_backup_ + B3_*m1_;

//      m1_backup_ = m1_;
//      m1_ = B3_.transpose()*m1_;
//      m0_ = B0_.transpose()*m0_ + B1_.transpose()*m1_backup_;
    }
    timer.stop();

    std::cout << "P*m=\n" << m0_ << "\n" << m1_ << std::endl;
//    std::cout << "P'*m=\n" << m0_ << "\n" << m1_ << std::endl;
    printf("for %zu iterations using Eigen 3x3 blocks, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());
  }

  {
    //using Eigen 6x6 matrices
    Eigen::Matrix<Real,6,6> P_;
    Eigen::Matrix<Real,6,1> m_;

    //copy to Eigen matrices
    copyMat6bToArray(P,P_.data());
    copyVec6bToArray(m0,m_.data());

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
      m_ = P_*m_;
//      m_ = P_.transpose()*m_;
    }
    timer.stop();

    std::cout << "P*m=\n" << m_ << std::endl;
//    std::cout << "P'*m=\n" << m_ << std::endl;
    printf("for %zu iterations using Eigen 6x6 matrices, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

  }
}

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
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
