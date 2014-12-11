#include <test.h>
#include <../mytest.h>

//int main()
int main(int argc, char *argv[]) //use this for console output
{
	//in test.h
	//test_common();

	//test_linalg3();
	//test_transform();
	//test_spatial();
	//test_matrix();

	//test_surface();
	//test_updateWheelContactGeom();
	//test_updateTrackContactGeom();

	//test_stateToHT();
	//test_qvelToQdot();

	//test_wheelJacobians();
	//test_trackJacobians();
	//test_forwardVelKin();
	//test_initTerrainContact();

	//test_subtreeInertias();
	//test_jointSpaceInertia();
	//test_jointSpaceBiasForce();
	//test_forwardDyn();

	test_simulate();


	//in mytest.h

	//test_convertToWmrModelODE();
	//test_simulate_ODE();
	//test_benchmark();

	
	std::cin.get();
}



