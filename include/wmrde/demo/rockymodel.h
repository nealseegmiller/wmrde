#ifndef _WMRSIM_ROCKYMODEL_H_
#define _WMRSIM_ROCKYMODEL_H_

#include <state.h>
#include <algebra/matrix.h>
#include <wheelgroundcontact.h>
#include <actuator.h>

//& to pass object by reference
void rocky(WmrModel& mdl, Real state[], Real qvel[]);

void rockyController(const WmrModel& mdl, const Real time, const Real state[], //inputs
		Real u[], Real qvel_cmd[]); //outputs

inline void rockyAct( const Real params[], const Real ucmd[], const Real u[], const Real interr[], //inputs
		Real f[], Real err[], Real* dfdu) {

	const int na = 8; //6 independently actuated wheels, 2 steer angles
	setMat(na,na,0.0,dfdu);
	for (int i=0; i<na; i++)
		PIact( params, ucmd[i], u[i], interr[i], //inputs
				f[i], err[i], dfdu[S2I(i,i,na)]); //outputs
}

void rockyConstraints( const WmrModel& mdl, const Real jd[], const Real jr[], //inputs
		Real c[], Real Jc[], Real f[], Real df_djd[], Real df_djr[]); //outputs

#endif