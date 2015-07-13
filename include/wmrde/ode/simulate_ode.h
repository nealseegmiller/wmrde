//simulateODE.h
//functions to convert a WmrModel object to a WmrModelODE object (Open Dynamics Engine) and simulate it

#ifndef _WMRSIM_SIMULATEODE_H_
#define _WMRSIM_SIMULATEODE_H_

#include <state.h>
#include <ode/WmrModelODE.h>
#include <collision.h>


//pass by reference

void convertToWmrModelODE( const WmrModel &mdl, WmrModelODE &mdl_ode);

void getHTODE( const WmrModel &mdl, const WmrModelODE &mdl_ode, HomogeneousTransform HT_world[], HomogeneousTransform HT_parent[] );
void getStateODE( const WmrModelODE &mdl_ode, Real state[] ); 
void setStateODE( const WmrModel &mdl, const Real state[], WmrModelODE &mdl_ode );

//joint space velocity
void getQvelODE( const WmrModelODE &mdl_ode, Real qvel[] );
void setQvelODE( const WmrModel &mdl, const Real qvel[], WmrModelODE &mdl_ode );

void stepODE( const Real time, const WmrModel& mdl, const SurfaceVector& surfaces, const Real dt, //inputs
	WmrModelODE& mdl_ode, HomogeneousTransform HT_parent[], WheelContactGeom contacts[]); //outputs

#endif
