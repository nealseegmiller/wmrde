#ifndef _WMRDE_WHEELGEOM_H_
#define _WMRDE_WHEELGEOM_H_

#include <wmrde/algebra/transform.h>

namespace wmrde
{

//TODO, move this?
struct Tangent //TODO, rename this?
{
  Vec3 point; //!< point of tangency
  Vec3 direction; //!< direction of tangent line.
};

class WheelGeom
{
 public:
  Real radius; //!< radius of the wheel (or driving sprocket for track)
  std::vector<Tangent> tangents; /*!< Discretization of wheel surface into potential contact points.
  Each Tangent specifies point and the direction point moves if joint rate > 0, in wheel coords
  */

  //TODO, implement a virtual clone() function
};

} //namespace

#endif
