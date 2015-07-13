#ifndef _WMRSIM_TALONMODEL_H_
#define _WMRSIM_TALONMODEL_H_

#include <state.h>
#include <algebra/matrix.h>
#include <wheelgroundcontact.h>
#include <actuator.h>

//& to pass object by reference
void talon(WmrModel& mdl, Real state[], Real qvel[]);

void talonController(const WmrModel& mdl, const Real time, const Real state[], //inputs
		Real u[], Real qvel_cmd[]); //outputs

inline void talonAct( const Real params[], const Real ucmd[], const Real u[], const Real interr[], //inputs
		Real f[], Real err[], Real* dfdu) {

	const int na = 2; //2 independently actuated sprockets
	setMat(na,na,0.0,dfdu);
	for (int i=0; i<na; i++)
		PIact( params, ucmd[i], u[i], interr[i], //inputs
				f[i], err[i], dfdu[S2I(i,i,na)]); //outputs
}

#endif