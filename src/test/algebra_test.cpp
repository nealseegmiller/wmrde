#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <wmrde/algebra/linalg3.h>
#include <wmrde/algebra/matrix.h>

using namespace wmrde;

// Declare a test
TEST(TestSuite, addTwoIntsPass)
{
  EXPECT_EQ(addTwoInts(2,2),4);
}

// Declare another test
TEST(TestSuite, addTwoIntsFail)
{
  EXPECT_EQ(addTwoInts(2,2),5);
}

TEST(TestSuite, multMatMat3)
{
  Mat3 A,B,C,C2;

  setMat3(1.25,2,3,4,5,6,7,8,9,A);
  setMat3(9,8,7,6,5,4,3,2,1,B);

  multMatMat3(A,B,C);
  printf("C = \n%s", Mat3ToString(C).c_str());

  Eigen::Matrix<Real,3,3> A_;
  Eigen::Matrix<Real,3,3> B_;
  Eigen::Matrix<Real,3,3> C_;

  //.data() returns a pointer to the data array of the matrix
  copyMat3ToArray(A, A_.data());
  copyMat3ToArray(B, B_.data());
  C_ = A_*B_;
  copyArrayToMat3(C_.data(), C2);

  printf("C2 = \n%s", Mat3ToString(C2).c_str());

  EXPECT_TRUE(Mat3Equal(C,C2));
}
// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
