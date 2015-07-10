

#ifndef _WMRSIM_TEST_ODE_H_
#define _WMRSIM_TEST_ODE_H_

#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
//#include <memory> //for unique_ptr, C++11 only

#include <Eigen/Dense>

#include <demo/models.h>
#include <demo/terrains.h>
#include <dynamics.h>
#include <ode/simulateODE.h>

void test_convertToWmrModelODE();
void test_simulate_ODE();
void test_benchmark();


#endif
