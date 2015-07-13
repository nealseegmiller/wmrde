//actuator.h
//functions for the calculation of actuator forces.
//required for dynamic simulation.

#ifndef _WMRDE_ACTUATOR_H_
#define _WMRDE_ACTUATOR_H_

#include <wmrde/algebra/matrix.h>

void PIact( const Real params[], const Real ucmd, const Real u, const Real interr, //inputs
	Real& f, Real& err, Real& dfdu); //outputs

#endif
