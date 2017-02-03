#include <benchmark/benchmark.h>

#include <iostream>
#include <wmrde/algebra/rotation.h>
#include <wmrde/algebra/spatial.h>

using namespace wmrde;

typedef Eigen::Matrix<Real,6,6> Mat6;
typedef Eigen::Matrix<Real,6,1> Vec6;

// ------------------------------ 3x3 matrix multiplication ------------------------------
static void BM_MatMatMult3x3(benchmark::State& state) {
  Mat3 A, B, C;
  A.setRandom();
  B.setRandom();
//  std::cout << "A = \n" << A << std::endl;
//  std::cout << "B = \n" << B << std::endl;
  while (state.KeepRunning())
  {
    C.noalias() = A*B;
  }
}
BENCHMARK(BM_MatMatMult3x3); //2 ns

static void BM_MatTMatMult3x3(benchmark::State& state) {
  Mat3 A, B, C;
  A.setRandom();
  B.setRandom();
  while (state.KeepRunning())
  {
    C.noalias() = A*B;
  }
}
BENCHMARK(BM_MatTMatMult3x3); //2 ns

//------------------------------ rotation.h functions ------------------------------
inline double rand01() { return rand()/static_cast<double>(RAND_MAX); }
inline Real randAngle() { return 2.0*M_PI*(rand01() - 0.5); }

static void BM_Rotx(benchmark::State& state) {
  Mat3 Rx;
  Real angle = randAngle();
  while (state.KeepRunning())
  {
    Rx = Rotx(angle);
  }
}
BENCHMARK(BM_Rotx); //2 ns

static void BM_EulerToRot(benchmark::State& state) {
  Mat3 R;
  Real roll = randAngle();
  Real pitch = randAngle();
  Real yaw = randAngle();
//  printf("roll = %f, pitch = %f, yaw = %f\n", roll, pitch, yaw);

  while (state.KeepRunning())
  {
    R = eulerToRot(roll,pitch,yaw);
  }
}
BENCHMARK(BM_EulerToRot); //100 ns, why so slow?

static void BM_RotxUsingAngleAxis(benchmark::State& state) {
  Mat3 Rx;
  Real angle = randAngle();
  while (state.KeepRunning())
  {
    Rx = RotxTest(angle);
  }
}
BENCHMARK(BM_RotxUsingAngleAxis); // 2 ns

static void BM_EulerToRotUsingAngleAxis(benchmark::State& state) {
  Mat3 R;
  Real roll = randAngle();
  Real pitch = randAngle();
  Real yaw = randAngle();
  while (state.KeepRunning())
  {
    R = eulerToRotTest(roll,pitch,yaw);
  }
}
BENCHMARK(BM_EulerToRotUsingAngleAxis); //2 ns

// ------------------------------ transform.h functions ------------------------------
inline Mat3 randomRotation() { return eulerToRot(randAngle(),randAngle(),randAngle()); }
inline HTransform randomHTransform() { return HTransform(randomRotation(), Vec3::Random()); }

static void BM_HTcomposeWith(benchmark::State& state) {
  HTransform T = randomHTransform();
  HTransform TT;
  while (state.KeepRunning())
  {
    TT = T.composeWith(T);
  }
}
BENCHMARK(BM_HTcomposeWith); //2 ns

static void BM_HTcomposeInvWith(benchmark::State& state) {
  HTransform T = randomHTransform();
  HTransform TT;
  while (state.KeepRunning())
  {
    TT = T.composeInvWith(T);
  }
}
BENCHMARK(BM_HTcomposeInvWith); //2 ns

static void BM_HTcomposeWithUsing4x4(benchmark::State& state) {
  Eigen::Matrix<Real,4,4> T, TT;
  T = randomHTransform().to4x4();
  while (state.KeepRunning())
  {
    TT.noalias() = T*T;
  }
}
BENCHMARK(BM_HTcomposeWithUsing4x4); //10 ns

// ------------------------------ spatial.h matrix-vector multiplication ------------------------------
inline Mat6b randomPlucker() { return HTToPlucker(randomHTransform()); }
inline Vec6b randomVec6b() { return Vec6b(Vec3::Random(), Vec3::Random()); }

static void BM_multPluckerVec6b(benchmark::State& state) {
  Mat6b P = randomPlucker();
  Vec6b v = randomVec6b();
  Vec6b w; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    w = multPluckerVec6b(P,v);
  }
}
BENCHMARK(BM_multPluckerVec6b); //2 ns

static void BM_multPluckerTVec6b(benchmark::State& state) {
  Mat6b P = randomPlucker();
  Vec6b v = randomVec6b();
  Vec6b w; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    w = multPluckerTVec6b(P,v);
  }
}
BENCHMARK(BM_multPluckerTVec6b); //2 ns

static void BM_multPluckerVecUsing6x6(benchmark::State& state) {
  Mat6 P = randomPlucker().to6x6();
  Vec6 v = Vec6::Random();
  Vec6 w; //result
  while (state.KeepRunning())
  {
    w.noalias() = P*v;
  }
}
BENCHMARK(BM_multPluckerVecUsing6x6); //2 ns

static void BM_multPluckerTVecUsing6x6(benchmark::State& state) {
  Mat6 P = randomPlucker().to6x6();
  Vec6 v = Vec6::Random();
  Vec6 w; //result
  while (state.KeepRunning())
  {
    w.noalias() = P.transpose()*v;
  }
}
BENCHMARK(BM_multPluckerTVecUsing6x6); //8 ns

// ------------------------------ spatial.h matrix-matrix multiplication ------------------------------
inline Mat6b randomInertia() { return toSpatialInertia(rand01(), Vec3::Random(), Mat3::Random()); }

static void BM_multPluckerTMat6b(benchmark::State& state) {
  Mat6b P = randomPlucker();
  Mat6b I = randomInertia();
  Mat6b M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M = multPluckerTMat6b(P,I);
  }
}
BENCHMARK(BM_multPluckerTMat6b); //49 ns

static void BM_multMat6bPlucker(benchmark::State& state) {
  Mat6b P = randomPlucker();
  Mat6b I = randomInertia();
  Mat6b M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M = multMat6bPlucker(I,P);
  }
}
BENCHMARK(BM_multMat6bPlucker); //65 ns

static void BM_multPluckerTMatUsing6x6(benchmark::State& state) {
  Mat6 P = randomPlucker().to6x6();
  Mat6 I = randomInertia().to6x6();
  Mat6 M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M.noalias() = P.transpose()*I;
  }
}
BENCHMARK(BM_multPluckerTMatUsing6x6); //51 ns

static void BM_multMatPluckerUsing6x6(benchmark::State& state) {
  Mat6 P = randomPlucker().to6x6();
  Mat6 I = randomInertia().to6x6();
  Mat6 M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M.noalias() = I*P;
  }
}
BENCHMARK(BM_multMatPluckerUsing6x6); //33 ns

// ------------------------------ dynamic vs. fixed matrix operations ------------------------------

template<int T, int N>
static void BM_MatOp(benchmark::State& state) {
  Eigen::Matrix <Real,T,T> A, B, C;
  A.resize(N,N);
  B.resize(N,N);
  C.resize(N,N);
  A.setRandom();
//  B.setRandom();
  B = A.transpose()*A; //symmetric

//  std::cout << "A = \n" << A << std::endl;
//  std::cout << "B = \n" << B << std::endl;

  while (state.KeepRunning())
  {
//    C.noalias() = A*B;
//    C.noalias() = A.transpose()*B;

    Eigen::LLT< Eigen::Matrix<Real,T,T> > A_llt(A);
    C = A_llt.matrixL();
  }
}
BENCHMARK_TEMPLATE(BM_MatOp, 12, 12); //24 ns
BENCHMARK_TEMPLATE(BM_MatOp, Eigen::Dynamic, 12); //226 ns

BENCHMARK_MAIN();
