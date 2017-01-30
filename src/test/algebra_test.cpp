#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <wmrde/algebra/linalg3.h>
#include <wmrde/util/timer.h>

using namespace wmrde;


TEST(TestSuite, addVec3)
{
  Vec3 a,b,c,c2;
  Real eps = 1e-6;
  setVec3(eps, a);
  setVec3(eps, b);

  Vec3 a_backup; copyVec3(a,a_backup);

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
    copyVec3(a_backup,a);

    Eigen::Matrix<Real,3,1> a_, b_, c_;
    copy3(a, a_.data());
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

  Mat3 A_backup; copyMat3(A,A_backup);

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
    copyMat3(A_backup,A);

    Eigen::Matrix<Real,3,3> A_, B_, C_;
    copyMat3ToArray(A, A_.data());
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

  Mat3 A_backup; copyMat3(A,A_backup);

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
    copyMat3(A_backup,A);

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

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
