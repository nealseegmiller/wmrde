#ifndef __WMRSIM_TALONSCENE_H_
#define __WMRSIM_TALONSCENE_H_

#include <WmrModel.h>
#include <animate/WmrAnimation.h>
#include <contactgeom.h>

void talonScene(const WmrModel& mdl, const TrackContactGeom contacts[], WmrAnimation& anim ); 

#endif