//functions to generate random transforms for unit testing and
//benchmarking

#ifndef _WMRDE_RANDOM_H_
#define _WMRDE_RANDOM_H_

#include <wmrde/algebra/random.h>
#include <wmrde/algebra/spatial.h>

namespace wmrde
{

//TODO, function to set seed of random number generator
inline double rand01() { return rand()/static_cast<double>(RAND_MAX); }
inline Real randAngle() { return 2.0*M_PI*(rand01() - 0.5); }
inline Mat3 randomRotation() { return eulerToRot(randAngle(),randAngle(),randAngle()); }
inline HTransform randomHTransform() { return HTransform(randomRotation(), Vec3::Random()); }
inline Mat6 randomPlucker() { return HTToPlucker(randomHTransform()); }
inline Mat6 randomInertia()
{
  Mat3 Rand3 = Mat3::Random();
  return toSpatialInertia(rand01(), Vec3::Random(), Rand3*Rand3.transpose());
}

} //namespace

#endif //_WMRDE_RANDOM_H_
