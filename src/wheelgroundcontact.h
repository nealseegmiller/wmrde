//wheelgroundcontact.h
//functions for the computation of wheel-ground contact forces
//required for dynamic simulation

#ifndef _WMRSIM_WHEELGROUNDCONTACT_H_
#define _WMRSIM_WHEELGROUNDCONTACT_H_

//wheel-ground contact models

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

#include <algebra/linalg3.h>
#include <common/interp.h>


void setWgcParams( const Real Kp, Real params[] );

void odeWgc( const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
		Vec3 fw, Real J[]); //outputs

void pacejkaWgc( const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
		Vec3 fw, Real J[]); //outputs

void ishigamiWgc(const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
    Vec3 fw, Real J[]); //outputs

void ishigamiLUTWgc( const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
		Vec3 fw, Real J[]); //outputs

void pacejkaDerivatives(const Real Ba, const Real Bs, const Real Ca, const Real Cs, const Real Da, const Real Ds, const Real Ea, const Real Es, const Real Ka, const Real Ks, const Real alpha, const Real fz, const Real mu, const Real s, const Real sgn, //inputs
			Real* dfx_ds, Real* dfy_ds, Real* dfx_dalpha, Real* dfy_dalpha); //outputs

void calcSlip( const Real vx, const Real vy, const Real Rw, const int method, //inputs
		Real* s, Real* alpha, Real* ds_dvx, Real* ds_dRw, Real* dalpha_dvx, Real* dalpha_dvy, Real* dalpha_dRw); //outputs

#define ODE_WGC 0
#define PACEJKA_WGC 1
#define ISHIGAMI_WGC 2
#define ISHIGAMI_LUT_WGC 3
//modify this:
#define WGC_MODEL_TYPE ISHIGAMI_LUT_WGC

//uniform wheel-ground contact model. assumes all wheels are identical (wheelno is not used)
inline void uniformWgc( const int wheelno, const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
				Vec3 fw, Real J[]) { //outputs
	
#if WGC_MODEL_TYPE == ODE_WGC
	odeWgc(params,vc,Rw,dz, fw,J);
#elif WGC_MODEL_TYPE == PACEJKA_WGC
	pacejkaWgc(params,vc,Rw,dz, fw,J);
#elif WGC_MODEL_TYPE == ISHIGAMI_WGC
	ishigamiWgc(params,vc,Rw,dz, fw,J);
#elif WGC_MODEL_TYPE == ISHIGAMI_LUT_WGC
	ishigamiLUTWgc(params,vc,Rw,dz, fw,J);
#endif
}


#endif
