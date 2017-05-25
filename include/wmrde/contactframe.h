#ifndef _WMRDE_CONTACTFRAME_H_
#define _WMRDE_CONTACTFRAME_H_

#include <wmrde/algebra/transform.h>

namespace wmrde
{

/*!
 * ContactFrame orientation convention:
 * z-axis is ground surface normal
 * x-axis is orthogonal to wheel frame y-axis
 */
struct ContactFrame
{
  int parent_idx; //!< index of parent wheel frame in WmrModel
  Vec3 wheel_tangent_dir; //!< tangent direction to wheel in wheel coords, the direction point moves if joint rate > 0
  HTransform HT_wheel; //!< transform from contact frame to wheel coords.
  HTransform HT_world; //!< transform from contact frame to world coords.

  Real dz; //!< contact height error. < 0 indicates contact point is below ground surface
  bool inContact() const { return dz < 0; }
};

} //namespace

#endif
