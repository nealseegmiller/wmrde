#ifndef _WMRDE_ZOEMODEL_H_
#define _WMRDE_ZOEMODEL_H_

#include <wmrde/wmrmodel.h>
#include <visualization_msgs/MarkerArray.h>

namespace wmrde
{

/*!
 * Initialize WmrModel to Zoe rover design.
 * \param mdl the WmrModel object to initialize
 * \param ma marker array containing mesh resource markers for visualization in rviz.
 */
void makeZoeModel(WmrModel& mdl);

//TODO, move this to its own file to separate out ros dependency?
void makeZoeModelMarkers(
    const WmrModel& mdl,
    visualization_msgs::MarkerArray& markers);

} //namespace

/*
//& to pass object by reference
void zoe(WmrModel& mdl, Real state[], Real qvel[]);

void zoeController(const WmrModel& mdl, const Real time, const Real state[], //inputs
		Real u[], Real qvel_cmd[]); //outputs

void zoeArcController( const Real speed, const Real turnrad, const Real* steer, //inputs
	Real* steer_cmd, Real yawrate, Real* wheelvel, Real* steer_rates); //outputs

inline void zoeAct( const Real params[], const Real ucmd[], const Real u[], const Real interr[], //inputs
		Real f[], Real err[], Real* dfdu) {

	const int na = 4; //4 independently actuated wheels
	setMat(na,na,0.0,dfdu);
	for (int i=0; i<na; i++)
		PIact( params, ucmd[i], u[i], interr[i], //inputs
				f[i], err[i], dfdu[S2I(i,i,na)]); //outputs
}

void zoeConstraints( const WmrModel& mdl, const Real jd[], const Real jr[], //inputs
		Real c[], Real Jc[], Real f[], Real df_djd[], Real df_djr[]); //outputs
*/

#endif
