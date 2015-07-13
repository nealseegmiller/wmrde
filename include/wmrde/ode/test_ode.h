

#ifndef _WMRDE_TEST_ODE_H_
#define _WMRDE_TEST_ODE_H_

#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
//#include <memory> //for unique_ptr, C++11 only

#include <Eigen/Dense>

#include <wmrde/demo/models.h>
#include <wmrde/demo/terrains.h>
#include <wmrde/dynamics.h>
#include <wmrde/ode/simulate_ode.h>

void test_convertToWmrModelODE();
void test_simulate_ODE();
void test_benchmark();


#endif
