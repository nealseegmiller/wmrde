//eigensolve.h
//functions for solving linear systems and decomposing matrices using the Eigen library.
//use fixed size matrices for faster computation.
//if not included, must include alternative files to define:
//solve(), subset(), and chol()

#ifndef _WMRDE_EIGENSOLVE_H_
#define _WMRDE_EIGENSOLVE_H_

#include <Eigen/Dense>
#include <wmrde/algebra/matrix.h>

//print matrix size if dynamic (vs. fixed)
#define PRINT_MATRIX_SIZE_IF_DYNAMIC 0

//Eigen is faster for fixed size matrices
//to used fixed size matrices for solve(), subset()

//zoe
//#define FIXED_NROWS_0 12+0 /*+0 if fix front axle roll, else +1*/
//#define FIXED_NCOLS_0 9+0

//rocky
//#define FIXED_NROWS_0 19
//#define FIXED_NCOLS_0 10

//#define FIXED_NROWS_1
//#define FIXED_NCOLS_1

//for chol()

//zoe
//#define FIXED_N_0 13+0
//#define FIXED_N_1 9+0 /*ideal actuators*/

//rocky
//#define FIXED_N_0 18
//#define FIXED_N_1 10 /*ideal actuators*/

typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;

//A*x = b, solve for x
void eigenSolveDynamic( const int nrows, const int ncols, Real* A, Real* b, Real* x );
void eigenSolveFixed( const int nrows, const int ncols, Real* A, Real* b, Real* x );
//call this one:
inline void solve(const int nrows, const int ncols, Real* A, Real* b, Real* x) {
	eigenSolveFixed(nrows,ncols,A,b,x);	
}

//compute the Cholesky decomposition of A = LL^T, L is lower triangular
//return true if Success. to succeed A must be positive definite
bool eigenCholDynamic( const int n, Real* A, Real* L);
bool eigenCholFixed( const int n, Real* A, Real* L);
//call this one:
inline bool chol( const int n, Real* A, Real* L) {
	return eigenCholFixed(n,A,L);
}

#endif  //_WMRDE_EIGENSOLVE_H_
