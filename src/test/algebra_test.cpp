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
  Vec3 a,b,c,c2;
  Real eps = 1e-6;
  setVec3(eps, a);
  setVec3(eps, b);

  Vec3 a0; copyVec3(a,a0);

  size_t num_iter = 1e7;
  Timer timer; timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    addVec3(a,b,c);
    copyVec3(c,a);
//    asm(""); //empty assembly to prevent gcc from optimizing out the loop
  }
  timer.stop();

  printf("for %zu iterations of addVec3, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

  //verify with Eigen library
  {
    Eigen::Matrix<Real,3,1> a_, b_, c_;
    copy3(a0, a_.data());
    copy3(b, b_.data());

    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      c_ = a_ + b_;
      a_ = c_;
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of Eigen 3x1 vector addition, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copy3(c_.data(), c2);
    printf("c = \n%s", Vec3ToString(c).c_str());
    printf("c2 = \n%s", Vec3ToString(c2).c_str());

    EXPECT_TRUE(Vec3Equal(c,c2));
  }
}

TEST(TestSuite, addMat3)
{
  Mat3 A,B,C,C2;
  Real eps = 1e-6;
  setMat3Diagonal(eps, eps, eps, A);
  setMat3Diagonal(eps, eps, eps, B);

  Mat3 A0; copyMat3(A,A0); //backup

  size_t num_iter = 1e7;
  Timer timer; timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    addMat3(A,B,C);
    copyMat3(C,A);
//    asm(""); //empty assembly to prevent gcc from optimizing out the loop
  }
  timer.stop();

  printf("for %zu iterations of addMat3, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

  //verify with Eigen library
  {
    Eigen::Matrix<Real,3,3> A_, B_, C_;
    copyMat3ToArray(A0, A_.data());
    copyMat3ToArray(B, B_.data());

    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      C_ = A_ + B_;
      A_ = C_;
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of Eigen 3x3 matrix addition, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copyArrayToMat3(C_.data(), C2);
    printf("C = \n%s", Mat3ToString(C).c_str());
    printf("C2 = \n%s", Mat3ToString(C2).c_str());

    EXPECT_TRUE(Mat3Equal(C,C2));
  }
}

TEST(TestSuite, multMatMat3)
{
  Mat3 A,B,C,C2;
  Real eps = 1e-6;
  setMat3Diagonal(1+eps, 1+eps, 1+eps, A);
  setMat3Diagonal(1+eps, 1+eps, 1+eps, B); B[1] = eps;

  Mat3 A0; copyMat3(A,A0);

  size_t num_iter = 1e7;
  Timer timer; timer.start();
  for (size_t iter = 0; iter < num_iter; iter++)
  {
    multMatMat3(A,B,C);
    copyMat3(C,A);
//    asm(""); //empty assembly to prevent gcc from optimizing out the loop
  }
  timer.stop();

  printf("for %zu iterations of multMatMat3, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

  //verify with Eigen library
  {
    copyMat3(A0,A);

    Eigen::Matrix<Real,3,3> A_, B_, C_;
    copyMat3ToArray(A, A_.data());
    copyMat3ToArray(B, B_.data());

    timer.start();
    for (size_t iter = 0; iter < num_iter; iter++)
    {
      C_ = A_*B_;
      A_ = C_;
//      asm(""); //empty assembly to prevent gcc from optimizing out the loop
    }
    timer.stop();
    printf("for %zu iterations of Eigen 3x3 matrix multiplication, elapsed time = %f ms\n", num_iter, timer.elapsedTimeMs());

    copyArrayToMat3(C_.data(), C2);
    printf("C = \n%s", Mat3ToString(C).c_str());
    printf("C2 = \n%s", Mat3ToString(C2).c_str());

    EXPECT_TRUE(Mat3Equal(C,C2));
  }
}

TEST(TestSuite, composeHT)
{
  Timer timer;

  VecEuler euler;
  Vec3 translation;
  HomogeneousTransform HT,HT2;

  setEuler(DEGTORAD(0),DEGTORAD(0),DEGTORAD(0.1),euler);
  Real eps = 1e-6;
  setVec3(eps,eps,eps,translation);
  poseToHT(euler,translation,HT);

  HomogeneousTransform HT0; copyHT(HT, HT0); //backup

  std::cout << "HT=\n"; printHT(HT,-1,-1);

  //test compose HT
  int num_iter = (int) 1e7;
  {
    //time it
    timer.start();
    for (int i=0; i<num_iter; i++)
    {
      composeHT(HT,HT0,HT2);
      copyHT(HT2,HT);
    }
    timer.stop();

    std::cout << "HT*HT=\n"; printHT(HT2,-1,-1);
    std::cout << "iterations: " << (Real) num_iter << std::endl;
    std::cout << "elapsed time (ms): " << timer.elapsedTimeMs() << std::endl;
  }

  //validate using Eigen::Matrix<Real,3,3>
  {
    Eigen::Matrix<Real,3,3> R0,R;
    Eigen::Matrix<Real,3,1> t0,t;

    copyMat3ToArray(HT0,R.data());
    copy3(HT0+COL3,t.data());
    std::cout << "Eigen HT=\n" << R << "\n" << t << std::endl;

    R0 = R; //backup
    t0 = t;

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
      t = R*t0 + t;
      R = R*R0;
    }
    timer.stop();

    std::cout << "Eigen HT*HT=\n" << R << "\n" << t << std::endl;
    std::cout << "iterations: " << (Real) num_iter << std::endl;
    std::cout << "elapsed time (ms): " << timer.elapsedTimeMs() << std::endl;
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

    std::cout << "Eigen 4x4 HT*HT=\n" << HT_ << std::endl;
    std::cout << "iterations: " << (Real) num_iter << std::endl;
    std::cout << "elapsed time (ms): " << timer.elapsedTimeMs() << std::endl;
  }

}

TEST(TestSuite, multPluckerVec6b)
{
  Timer timer;

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

  //inverse Plucker transform
  Mat6b P2;

  Vec6b m,f,v,v2;
  setVec6b(1,2,3,4,5,6,m);
//  setVec6b(2,3,4,5,6,7,f);

  std::cout << "m=\n"; printVec6b(m,-1,-1);
//  std::cout << "f=\n"; printVec6b(f,-1,-1);

  Vec6b m0; copyVec6b(m,m0);

  int num_iter = 1e7;

  //multiply vector by Plucker transform
  {
    //time it

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
//      multPluckerVec6b(P,m,v);
//      copyVec6b(v,m);
      multPluckerTVec6b(P,m,v2);
      copyVec6b(v2,m);
    }
    timer.stop();

    std::cout << "using spatial.h\n";
//    std::cout << "P*m=\n"; printVec6b(v,-1,-1);
    std::cout << "P'*m=\n"; printVec6b(v2,-1,-1);
    std::cout << "iterations: " << (Real) num_iter << std::endl;
    std::cout << "clock (ms): " << timer.elapsedTimeMs() << std::endl;
    std::cout << std::endl;
  }

  {
    //using Eigen
    Eigen::Matrix<Real,3,3> B0_, B1_, B3_;
    Eigen::Matrix<Real,3,1> m0_, m1_, m0_backup_, m1_backup_;

    //copy to Eigen matrices
    copyMat3ToArray(P,B0_.data());
    copyMat3ToArray(P+BLOCK1, B1_.data());
    copyMat3ToArray(P+BLOCK3, B3_.data());
    copy3(m0,m0_.data());
    copy3(m0+VEC3_SIZE, m1_.data());

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
//      m0_backup_ = m0_;
//      m0_ = B0_*m0_;
//      m1_ = B1_*m0_backup_ + B3_*m1_;

      m1_backup_ = m1_;
      m1_ = B3_.transpose()*m1_;
      m0_ = B0_.transpose()*m0_ + B1_.transpose()*m1_backup_;
    }
    timer.stop();

    std::cout << "using Eigen 3x3 blocks\n";
//    std::cout << "P*m=\n" << m0_ << "\n" << m1_ << std::endl;
    std::cout << "P'*m=\n" << m0_ << "\n" << m1_ << std::endl;
    std::cout << "iterations: " << (Real) num_iter << std::endl;
    std::cout << "clock (ms): " << timer.elapsedTimeMs() << std::endl;
  }

  {
    //using Eigen
    Eigen::Matrix<Real,6,6> P_;
    Eigen::Matrix<Real,6,1> m_;

    //copy to Eigen matrices
    copyMat6bToArray(P,P_.data());
    copyVec6bToArray(m0,m_.data());

    timer.start();
    for (int i=0; i<num_iter; i++)
    {
//      m_ = P_*m_;
      m_ = P_.transpose()*m_;
    }
    timer.stop();

    std::cout << "using Eigen\n";
//    std::cout << "P*m=\n" << m_ << std::endl;
    std::cout << "P'*m=\n" << m_ << std::endl;
    std::cout << "iterations: " << (Real) num_iter << std::endl;
    std::cout << "clock (ms): " << timer.elapsedTimeMs() << std::endl;
  }
}

TEST(TestSuite, matrix)
{
  //matrix.h vs. Eigen

  //for timing
  int n = (int) 1e6;
  Timer timer;

  const int SIZE = 16;

  //matrix.h matrices
  Real A[SIZE*SIZE];
  Real B[SIZE*SIZE];
  Real B0[SIZE*SIZE]; //backup, to initialize Eigen
  Real B2[SIZE*SIZE];
  Real b[SIZE];

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


  //matrix.h
  int ri = 2;
  int ci = 3;
  int brows = SIZE/2;
  int bcols = SIZE/2;
  Real val = 2.0;
  Real eps = 1e-6;
//  Real m = 1.0 + eps;
  Real m = 1.0;

  timer.start();
  for (int i=0; i<n; i++)
  {
//    setMat(SIZE,SIZE,val,B);
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

//    addmMat(SIZE,SIZE,A,m,B);
    //addmMatRow(SIZE,SIZE,ri,A,SIZE,ri,m,B);
    //addmMatCol(SIZE,ci,A,ci,m,B);
//    addmMatBlock(SIZE,ri,ci,brows,bcols,A, SIZE,ri,ci,m,B);

//    multMatVec(SIZE,SIZE,A,b,m,B);
//    multMatTVec(SIZE,SIZE,A,b,m,B);
    multMatMat(SIZE,SIZE,A,SIZE,A,m,B);
//    multMatTMat(SIZE,SIZE,A,SIZE,A,m,B);
  }
  timer.stop();

//  printf("B[0] = %f\n", B[0]); //DEBUGGING
  std::cout << "(matrix.h) B=\n"; printMatReal(SIZE,SIZE,B,0,-1);
  std::cout << "iterations: " << (Real) n << std::endl;
  std::cout << "elapsed_time (ms): " << timer.elapsedTimeMs() << std::endl;
  std::cout << std::endl;

  if (1)
  {
    //Eigen

    //dynamic Eigen matrices
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> A_(SIZE,SIZE);
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> B_(SIZE,SIZE);
    Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> b_(SIZE,1);

    //fixed Eigen matrices
//    Eigen::Matrix<Real,SIZE,SIZE> A_;
//    Eigen::Matrix<Real,SIZE,SIZE> B_;
//    Eigen::Matrix<Real,SIZE,1> b_;

    memcpy(A_.data(), A, sizeof(Real)*SIZE*SIZE);
    memcpy(B_.data(), B0, sizeof(Real)*SIZE*SIZE);
    memcpy(b_.data(), b, sizeof(Real)*SIZE);

    //std::cout << "Eigen:\n";
    //std::cout << "A_=\n" << A_ << std::endl;
    //std::cout << "B_=\n" << B_ << std::endl;
    //std::cout << "b_=\n" << b_ << std::endl;

    timer.start();
    for (int i=0; i<n; i++)
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

//      B_.col(0) = A_*b_;
//      B_.col(0) = A_.transpose()*b_;
      B_ = A_*A_;
//      B_ = A_.transpose()*A_;
    }
    timer.stop();

    std::cout << "(Eigen) B=\n" << B_ << std::endl;
    std::cout << "iterations: " << (Real) n << std::endl;
    std::cout << "elapsed_time (ms): " << timer.elapsedTimeMs() << std::endl;
  }
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
