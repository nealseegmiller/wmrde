#include <gtest/gtest.h>

#include <wmrde/util/gridded_interpolant.h>
#include <wmrde/util/string_format.h>
#include <Eigen/Geometry> //to print the values of X,V

//uncomment one to enable/disable print output
#define DEBUG_INFO(x) do { std::cout << x << std::endl; } while(0)
//#define DEBUG_INFO(x) do { } while(0)

using namespace wmrde;

Eigen::Matrix<Real,-1,-1> GridVectorToMat(const GridVector& X)
{
  Eigen::Matrix<Real,-1,-1> M(X.numel,1);
  for (int i = 0; i < X.numel; i++)
  {
    M(i,0) = X.at(i);
  }
  return M;
}
Eigen::Matrix<Real,-1,-1> ArrayToMat(const Real* V,
    const int n0, //number of elements in first dimension
    const int n1=1) // " in 2nd dim
{
  Eigen::Matrix<Real,-1,-1> M(n0,n1);
  for (int i = 0; i < n0; i++)
  {
    for (int j = 0; j < n1; j++)
    {
      M(i,j) = V[i + j*n0];
    }
  }
  return M;
}

inline Real nan() { return std::numeric_limits<Real>::quiet_NaN(); }

TEST(TestSuite, linear_interp)
{
  const int N = 1;
  const int n0 = 5;
  std::array<GridVector,N> X = { GridVector(1.0, 1.0, n0) };
  std::vector<Real> V(n0);

  for (int i = 0; i < n0; i++)
  {
    V[i] = 2.0*i;
  }

  DEBUG_INFO("X[0] = \n" << GridVectorToMat(X[0]));
  DEBUG_INFO("V = \n" << ArrayToMat(V.data(),n0));

  GriddedInterpolant<Real,N> G(X,V,nan());

  std::array<Real,N> xq = {2.5};
  EXPECT_EQ(3.0, G.interpolate(xq));
}


TEST(TestSuite, bilinear_interp)
{
  const int N = 2;
  const int n0 = 5;
  const int n1 = 5;
  std::array<GridVector,N> X = {
      GridVector(1.0, 1.0, n0),
      GridVector(2.0, 2.0, n1)};

  GriddedInterpolant<Real,N> G(X, std::vector<Real>() ,nan());

  std::vector<Real> V(n0 * n1);
  for (int i = 0; i < n0; i++)
  {
    for (int j = 0; j < n1; j++)
    {
      G.getV(G.sub2ind(i,j)) = 2.0*i + 3.0*j;
    }
  }

  DEBUG_INFO("X[0] = \n" << GridVectorToMat(X[0]));
  DEBUG_INFO("X[1] = \n" << GridVectorToMat(X[1]));
  DEBUG_INFO("V = \n" << ArrayToMat(&G.getV(0),n0,n1));

  std::array<Real,N> xq = {2.5, 5.0};
  EXPECT_EQ(7.5,G.interpolate(xq));
}

TEST(TestSuite, trilinear_interp)
{
  const int N = 3;
  const int n0 = 5;
  const int n1 = 5;
  const int n2 = 3;
  std::array<GridVector,N> X = {
      GridVector(1.0, 1.0, n0),
      GridVector(2.0, 2.0, n1),
      GridVector(3.0, 3.0, n2)};

  GriddedInterpolant<Real,N> G(X, std::vector<Real>() ,nan());

  std::vector<Real> V(n0 * n1 * n2);
  for (int i = 0; i < n0; i++)
  {
    for (int j = 0; j < n1; j++)
    {
      for (int k = 0; k < n2; k++)
      {
        G.getV(G.sub2ind(i,j,k)) = 2.0*i + 3.0*j + 4.0*k;
      }
    }
  }

  //DEBUGGING
  for (int dim = 0; dim < N; dim++)
  {
    DEBUG_INFO("X[" << dim << "] = \n" << GridVectorToMat(X[dim]));
  }
  for (int k = 0; k < n2; k++)
  {
    Real* V_ptr = &G.getV(G.sub2ind(0,0,k));
    DEBUG_INFO("V[k=" << k << "] = \n" << ArrayToMat(V_ptr,n0,n1));
  }

  std::array<Real,N> xq = {2.5, 5.0, 4.5};
  EXPECT_EQ(9.5,G.interpolate(xq));
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
