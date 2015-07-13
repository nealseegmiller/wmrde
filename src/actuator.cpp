#include <wmrde/actuator.h>

void PIact( const Real params[], const Real ucmd, const Real u, const Real interr, //inputs
		   Real& f, Real& err, Real& dfdu) { //outputs
	Real Kp = params[0];
	Real Ki = params[1];

	Real fmax = params[2];

	err = ucmd-u;
	f = Kp*err + Ki*interr;
	if (fabs(f) > fmax) f=REALSIGN(f)*fmax; //enforce limit

	dfdu = -Kp;
	if (fabs(f) == fmax) dfdu=0.0;

}
