//actuator.h
//functions for the calculation of actuator forces.
//required for dynamic simulation.

#ifndef _WMRSIM_ACTUATOR_H_
#define _WMRSIM_ACTUATOR_H_

#include <algebra/matrix.h>

void PIact( const Real params[], const Real ucmd, const Real u, const Real interr, //inputs
	Real& f, Real& err, Real& dfdu); //outputs

#endif