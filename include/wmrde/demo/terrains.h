//terrains.h
//functions for the construction of various Surface objects
//all functions emplace a new surface at the back of a surface vector

#ifndef _WMRSIM_TERRAINS_H_
#define _WMRSIM_TERRAINS_H_

#include <vector>
#include <surface/PlaneSurf.h>
#include <surface/TriMeshSurf.h>
#include <surface/GridSurf.h>

void flat(SurfaceVector& surfaces);
void ramp(SurfaceVector& surfaces);
void grid(SurfaceVector& surfaces, const std::string FileName);


#endif