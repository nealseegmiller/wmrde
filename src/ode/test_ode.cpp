#include <ode/test_ode.h>

inline double tosec(timeval tim)
{
  return tim.tv_sec + (tim.tv_usec/1000000.0);
}

void test_convertToWmrModelODE() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	Real qvel[MAXNV];

	zoe(mdl,state,qvel);
	//rocky(mdl,state,qvel);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int ns = NUMSTATE(nf);
	const int nv = NUMQVEL(nf);

	setEuler(DEGTORAD(0), DEGTORAD(10), DEGTORAD(0), state+SI_ORIENT); //DEBUGGING

	//convert to WmrModelODE object
	WmrModelODE mdl_ode;
	convertToWmrModelODE(mdl, mdl_ode);

	setStateODE(mdl,state,mdl_ode);
	setQvelODE(mdl,qvel,mdl_ode);
	
	//test getStateODE
	std::cout << "state=\n"; printMatReal(ns,1,state,-1,-1);

	Real state_[MAXNS];
	getStateODE(mdl_ode,state_);
	std::cout << "getStateODE()=\n"; printMatReal(ns,1,state_,-1,-1); std::cout << std::endl;

	//test getQvelODE
	std::cout << "qvel=\n"; printMatReal(nv,1,qvel,-1,-1);

	Real qvel_[MAXNV];
	getQvelODE(mdl_ode,qvel_);
	std::cout << "getQvelODE()=\n"; printMatReal(nv,1,qvel_,-1,-1); std::cout << std::endl;

	if (1) {
		//DEBUGGING, center of mass offset
		int i = 0;

		Vec3 pos;
		mdl_ode.get_pos(i,pos);
		std::cout << "get_pos()=\n"; printVec3(pos,-1,-1);

		const dBodyID* body = mdl_ode.get_body();
		const dReal* pos_ = dBodyGetPosition(body[i]);
		std::cout << "dBodyGetPosition()=\n"; printVec3(pos_,-1,-1); std::cout << std::endl;
	}

	std::cout << std::endl;

}


void test_simulate_ODE() {

	//options
	bool do_anim = false;

	Real dt = .04;
	int nsteps = (int) floor(10.0/dt);
	Real time = 0;

	//constants
	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	Real qvel[MAXNV]; //for dynamic sim

	zoe(mdl,state,qvel);
	//rocky(mdl,state,qvel);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nw = mdl.get_nw();
	const int ns = NUMSTATE(nf);
	const int nv = NUMQVEL(nf);

	//terrain
	SurfaceVector surfs;

	//flat(surfs);
	//ramp(surfs);
	grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//contact geom
	WheelContactGeom wcontacts[WmrModel::MAXNW];

	initTerrainContact(mdl, surfs, wcontacts, state); //DEBUGGING

	//backup
	Real state0[MAXNS];
	copyVec(ns,state,state0); 

	//convert to WmrModelODE object
	WmrModelODE mdl_ode;
	convertToWmrModelODE(mdl, mdl_ode);

	setStateODE(mdl, state, mdl_ode);
	setQvelODE(mdl, qvel, mdl_ode);

	//allocate additional outputs
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	

#if WMRSIM_ENABLE_ANIMATION
	WmrAnimation anim;
	if (do_anim) { //animate
		anim.start();

		zoeScene(mdl, anim);
		//rockyScene(mdl, anim);

		for (int i=0; i<surfs.size(); i++)
			anim.addEntitySurface(surfs[i].get());
	}
#endif

	std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,state,-1,-1); std::cout << std::endl;

	for (int i=0; i<nsteps; i++) {

		stepODE( time, mdl, surfs, dt, mdl_ode, HT_parent, wcontacts);
		time += dt;

#if WMRSIM_ENABLE_ANIMATION
		if (do_anim) {
			anim.updateNodesLines(nf, HT_parent, nw, wcontacts);

			if (!anim.updateRender())
				goto stop;

			while (anim.get_mPause()) {
				if (!anim.updateRender())
					goto stop;
			}
		}
#endif
	}

	getStateODE(mdl_ode, state);
	std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,state,-1,-1); std::cout << std::endl;

#if WMRSIM_ENABLE_ANIMATION
	if (do_anim) {
		//DEBUGGING, wait for escape key before exiting
		while (1) {
			if (!anim.updateRender())
				goto stop;
		}
	}
	stop: //goto
#endif


	if (1) {
		//time it
		int n= (int) 100;
		clock_t t;

		time = 0;

//		t=clock();
    timeval t0,t1;
    gettimeofday(&t0, NULL);

		for (int iter=0; iter<n; iter++) {
			setStateODE(mdl, state0, mdl_ode);
			setQvelODE(mdl, qvel, mdl_ode);

			for (int i=0; i<nsteps; i++) {
				stepODE( time, mdl, surfs, dt, mdl_ode, HT_parent, wcontacts);
				time += dt;
			}
		}
		gettimeofday(&t1, NULL);
//		t=clock()-t;

		//getStateODE(mdl_ode, state);
		//std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,state,-1,-1); std::cout << std::endl;
		std::cout << "simulate ODE\n";
		std::cout << "iterations: " << (Real) n << std::endl;
//		std::cout << "clock (ms): " << t << std::endl;
		std::cout << "clock (ms): " << tosec(t1)-tosec(t0) << std::endl;
	}

}




void test_benchmark() {

	//options
	bool do_dyn = false;
	bool do_ode = false; //do Open Dynamics Engine

	const int nss = 4; //number of step sizes
	Real stepsizes[nss] = {.01, .02, .04, .1};

	//constants
	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	const int MAXNY = MAXNS + MAXNV + WmrModel::MAXNA; //for dynamic sim

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	Real qvel[MAXNV]; //for dynamic sim

	zoe(mdl,state,qvel);
	//rocky(mdl,state,qvel);

	//initialize wheel-ground contact model
	mdl.wheelGroundContactModel(0, mdl.wgc_p, 0, 0, 0, //inputs
					0, 0); //outputs
	
	mdl.wheelGroundContactModel=0; //use erp,cfm
	mdl.actuatorModel=0; //ideal actuators (using ideal actuators in Open Dynamics Engine)

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nw = mdl.get_nw();
	const int ns = NUMSTATE(nf); //number of states
	//for dynamic sim
	const int nv = NUMQVEL(nf); //size of joint space vel
	const int na = mdl.get_na(); 
	int ny;
	
	//convert to WmrModelODE object
	WmrModelODE mdl_ode;
	if (do_ode) {
		convertToWmrModelODE(mdl, mdl_ode);
	}

	//contact geom
	WheelContactGeom wcontacts[WmrModel::MAXNW];
	ContactGeom* contacts = static_cast<ContactGeom*>(wcontacts);
	
	//backup
	Real state0[MAXNS];
	copyVec(ns,state,state0); 

	//allocate
	Real y[MAXNY];
	Real ydot[MAXNY];

	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	
	clock_t tt[nss] = {0,0,0,0}; //total time, for each step size

	for (int surfno=0; surfno<10; surfno++) { //loop over surfaces

		std::cout << "surface number=" << surfno << std::endl;

		//vector of pointers to base class
		SurfaceVector surfs;
		std::string filename = ResourceDir() + std::string("gridsurfdata") + std::to_string(surfno+1) + std::string(".txt");
		grid(surfs, filename );

		copyVec(ns,state0,state); //DEBUGGING

		initTerrainContact(mdl, surfs, contacts, state);

		//init y
		if (do_dyn) {
			copyVec(ns,state,y);
			copyVec(nv,qvel,y+ns);
			setVec(na,0.0,y+ns+nv);
			ny = ns + nv + na;
		} else {
			copyVec(ns,state,y);
			ny = ns;
		}

		//backup
		Real y0[MAXNY];
		copyVec(ny,y,y0); 

		Real time = 0;
		std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,y,-1,-1); std::cout << std::endl;

		for (int ssno=0; ssno<nss; ssno++) { //loop over step sizes
	
			Real dt = stepsizes[ssno];
			const int nsteps = (int) floor(10.0/dt);

			std::cout << "dt=" << dt << std::endl;

			//time it
			int n= (int) 10;
			clock_t t;

			t=clock();
			for (int iter=0; iter<n; iter++) {

				time = 0;

				if (do_ode) {
					setStateODE(mdl, y0, mdl_ode);
					setQvelODE(mdl, qvel, mdl_ode);

					for (int i=0; i<nsteps; i++) {
						stepODE( time, mdl, surfs, dt, mdl_ode, HT_parent, wcontacts);
						time += dt;
					}
					getStateODE(mdl_ode,y);
				} else {
					copyVec(ny,y0,y);
					for (int i=0; i<nsteps; i++) {
						if (do_dyn) {
							odeDyn(time, y, mdl, surfs, contacts, dt, ydot, HT_parent);
						} else {
							odeKin(time, y, mdl, surfs, contacts, ydot, HT_parent);
						}
						addmVec(ny,ydot,dt,y); //Euler integration
						time += dt;
					}
				}
			}
			t=clock()-t;
			
			std::cout << "clock (ms): " << t << std::endl;

			std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,y,-1,-1); std::cout << std::endl;

			tt[ssno] += t; //total time
		}

	}

	for (int ssno=0; ssno<nss; ssno++) {
		std::cout << "step size: " << stepsizes[ssno] << ", clock (ms): " << tt[ssno] << std::endl;
	}
}
