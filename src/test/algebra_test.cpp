#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <wmrde/algebra/linalg3.h>
#include <wmrde/algebra/transform.h>
#include <wmrde/algebra/spatial.h>
#include <wmrde/algebra/matrix.h>
#include <wmrde/util/timer.h>

using namespace wmrde;


TEST(TestSuite, addVec3)
{
  printf("benchmark addVec3\n");

  Vec3 a,b,r,r2;
  Real eps = 1e-6;
  setVec3(eps, a);
  setVec3(eps, b);

  Vec3 a0; copyVec3(a,a0); //backup

  size_t num_iter = 1e7;
  Timer timer;
  {
    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      addVec3(a,b,a);
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of addVec3, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copyVec3(a,r); //first result
  }

  //verify with Eigen library
  {
    Eigen::Matrix<Real,3,1> a_, b_;
    copy3(a0, a_.data());
    copy3(b, b_.data());

    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      a_ = a_ + b_;
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of Eigen 3x1 vector addition, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copy3(a_.data(), r2); //second result
  }

  printf("r = \n%s", Vec3ToString(r).c_str());
  printf("r2 = \n%s", Vec3ToString(r2).c_str());

  EXPECT_TRUE(Vec3Equal(r,r2));
}

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

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
