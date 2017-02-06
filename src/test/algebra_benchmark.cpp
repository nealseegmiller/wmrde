#include <benchmark/benchmark.h>

#include <iostream>
#include <wmrde/algebra/rotation.h>
#include <wmrde/algebra/spatial.h>
#include <wmrde/algebra/spatial_block.h>
#include <wmrde/algebra/random.h>

using namespace wmrde;

//------------------------------ rotation.h functions ------------------------------
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
BENCHMARK(BM_EulerToRot); //2 ns

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

// ------------------------------ spatial.h convert HT to Plucker ------------------------------

static void BM_HTToPlucker(benchmark::State& state) {
  HTransform HT = randomHTransform();
  Mat6 P;

  while (state.KeepRunning())
  {
    P = HTToPlucker(HT);
  }
}
BENCHMARK(BM_HTToPlucker); // 2 ns

static void BM_invHTToPlucker(benchmark::State& state) {
  HTransform HT = randomHTransform();
  Mat6 P;

  while (state.KeepRunning())
  {
    P = invHTToPlucker(HT);
  }
}
BENCHMARK(BM_invHTToPlucker); // 2 ns

// ------------------------------ spatial.h Plucker-vector multiplication ------------------------------
static void BM_multPluckerVec(benchmark::State& state) {
  Mat6 P = randomPlucker();
  Vec6 v = Vec6::Random();
  Vec6 w; //result
  while (state.KeepRunning())
  {
    w.noalias() = multPluckerVec(P,v); //2 ns
  }
}
BENCHMARK(BM_multPluckerVec);

static void BM_multPluckerTVec(benchmark::State& state) {
  Mat6 P = randomPlucker();
  Vec6 v = Vec6::Random();
  Vec6 w; //result
  while (state.KeepRunning())
  {
    w.noalias() = multPluckerTVec(P,v);
  }
}
BENCHMARK(BM_multPluckerTVec); //2 ns

static void BM_multPluckerVecNonBlock(benchmark::State& state) {
  Mat6 P = randomPlucker();
  Vec6 v = Vec6::Random();
  Vec6 w; //result
  while (state.KeepRunning())
  {
    w.noalias() = P*v;
  }
}
BENCHMARK(BM_multPluckerVecNonBlock); //2 ns

static void BM_multPluckerTVecNonBlock(benchmark::State& state) {
  Mat6 P = randomPlucker();
  Vec6 v = Vec6::Random();
  Vec6 w; //result
  while (state.KeepRunning())
  {
    w.noalias() = P.transpose()*v;
  }
}
BENCHMARK(BM_multPluckerTVecNonBlock); //8 ns

// ------------------------------ spatial.h Plucker-inertia multiplication ------------------------------
static void BM_multPluckerTInertiaPlucker(benchmark::State& state) {
  Mat6 P = randomPlucker();
  Mat6 I = randomInertia();
  Mat6 M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M.noalias() = multPluckerTInertiaPlucker(P,I);
  }
}
BENCHMARK(BM_multPluckerTInertiaPlucker); //81 ns

static void BM_multPluckerTInertiaPluckerEigen(benchmark::State& state) {
  Mat6 P = randomPlucker();
  Mat6 I = randomInertia();
  Mat6 M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M.noalias() = P.transpose()*I*P;
  }
}
BENCHMARK(BM_multPluckerTInertiaPluckerEigen); //85 ns

//benchmark using deprecated Mat6b, requires spatial_block.h
static void BM_multPluckerTInertiaPlucker6b(benchmark::State& state) {
  Mat6b P; P.from6x6(randomPlucker());
  Mat6b I; I.from6x6(randomInertia());
  Mat6b M; //result

  HTransform TT;
  while (state.KeepRunning())
  {
    M = multMatPlucker6b(I,P);
    M = multPluckerTMat6b(P,M);
  }
}
BENCHMARK(BM_multPluckerTInertiaPlucker6b); //164 ns


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
    C.noalias() = A*B;
//    C.noalias() = A.transpose()*B;

//    Eigen::LLT< Eigen::Matrix<Real,T,T> > A_llt(A);
//    C = A_llt.matrixL();
  }
}
//time for matrix multiplication, time for cholesky decomposition
//BENCHMARK_TEMPLATE(BM_MatOp, 8, 8); //199 ns,
BENCHMARK_TEMPLATE(BM_MatOp, 12, 12); //500 ns, 24 ns
//BENCHMARK_TEMPLATE(BM_MatOp, 16, 16); //1014 ns,
//BENCHMARK_TEMPLATE(BM_MatOp, Eigen::Dynamic, 8); //274 ns,
BENCHMARK_TEMPLATE(BM_MatOp, Eigen::Dynamic, 12); //583 ns, 226 ns
//BENCHMARK_TEMPLATE(BM_MatOp, Eigen::Dynamic, 16); //1144 ns


// ------------------------------ main ------------------------------
BENCHMARK_MAIN();


