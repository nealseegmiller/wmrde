//terrains.h
//functions for the construction of various Surface objects
//all functions emplace a new surface at the back of a surface vector

#ifndef _WMRDE_TERRAINS_H_
#define _WMRDE_TERRAINS_H_

#include <vector>
#include <wmrde/surface/PlaneSurf.h>
#include <wmrde/surface/TriMeshSurf.h>
#include <wmrde/surface/GridSurf.h>

void flat(SurfaceVector& surfaces);
void ramp(SurfaceVector& surfaces);
void grid(SurfaceVector& surfaces, const std::string FileName);


#endif
