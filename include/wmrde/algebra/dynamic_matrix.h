//utility for linear algebra in R3

#ifndef _WMRDE_DYNAMIC_MATRIX_H_
#define _WMRDE_DYNAMIC_MATRIX_H_

#include <wmrde/util/real.h>
#include <Eigen/Dense>

namespace wmrde
{

typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> Matd;
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> Vecd;

} //namespace

#endif //_WMRDE_DYNAMIC_MATRIX_H_
