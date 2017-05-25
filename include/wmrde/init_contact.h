#ifndef _WMRDE_INIT_CONTACT_H_
#define _WMRDE_INIT_CONTACT_H_

#include <wmrde/wmrstate.h>
#include <wmrde/surface/surface.h>
#include <wmrde/collision.h>

namespace wmrde
{

inline void calcWmrModelContactFrames(
    const WmrModel& mdl,
    const Surfaces& surfs,
    const std::vector<HTransform>& HT_world,
    std::vector<ContactFrame>& contacts)
{
//  contacts.clear();
  contacts.resize(mdl.numWheelFrames());
  for (int i = 0; i < mdl.numWheelFrames(); i++)
  {
    int frame_idx = mdl.wheelFrameIndices()[i];
    calcContactFrameRootfind(surfs,
        *mdl.getFrame(frame_idx).wheel_geom, //TODO, check for nullptr
        HT_world[frame_idx],
        contacts[i]);
    contacts[i].parent_idx = frame_idx;
  }
}

//TODO, move this?
struct OptimizationOptions
{
  int max_iter;
  Real cost_tol;
  Real dcost_tol;

  OptimizationOptions() :
    max_iter(15),
    cost_tol(1e-6),
    dcost_tol(1e-6)
  {}
};

/*!
 * initialize contact of WmrModel with ground surface. The state is perturbed such
 * that the sum of squares of contact height error for all wheels is minimized.
 * \param[in] mdl the WmrModel
 * \param[in] surfs ground surfaces
 * \param[in,out] state the WmrState. the x, y, and yaw of body frame are held fixed as are actuated joint displacements
 * \param[out] contacts ContactGeom for every wheel
 * \param[in] options optional input to specify optimization accuracy
 */
void initWmrModelSurfaceContact(
    const WmrModel& mdl,
    const Surfaces& surfs,
    WmrState& state,
    std::vector<ContactFrame>& contacts,
    const OptimizationOptions& options = OptimizationOptions());

} //namespace

#endif
