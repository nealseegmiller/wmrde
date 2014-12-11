//eigensolve.h
//functions for solving linear systems and decomposing matrices using the Eigen library.
//use fixed size matrices for faster computation.
//if not included, must include alternative files to define:
//solve(), subset(), and chol()

#ifndef _WMRSIM_EIGENSOLVE_H_
#define _WMRSIM_EIGENSOLVE_H_

#include <Eigen/Dense>
#include <algebra/matrix.h>

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
//#define FIXED_N_1 16+0 /*use erp cfm*/
//#define FIXED_N_2 4 /*initTerrainContact*/

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

//obtain linearly independent subset of *rows* of A
//if is_ind[i] is true, row i belongs to the independent subset
void eigenSubsetDynamic( const int nrows, const int ncols, Real* A, const Real tol, bool is_ind[]);
void eigenSubsetFixed( const int nrows, const int ncols, Real* A, const Real tol, bool is_ind[]);
//call this one:
inline void subset( const int nrows, const int ncols, Real* A, const Real tol, bool is_ind[]) {
	eigenSubsetFixed(nrows,ncols,A,tol,is_ind);
}

//TODO, slow to use subfunction?
//templatize to avoid temporaries
//http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
template <typename DerivedA, typename DerivedB>
void eigenSetIsInd( const Eigen::MatrixBase<DerivedA>& R, const Eigen::MatrixBase<DerivedB>& I, const Real tol, bool is_ind[]) {
	//set is_ind bool array
	int Rrows = (int) R.rows();
	int Rcols = (int) R.cols();
	
	setVec(Rcols,false,is_ind); //initialize
	Real tol_ = tol*fabs(R(0,0)); //fabs is necessary!
	int n = std::min(Rrows,Rcols);
	for (int i=0; i<n; i++) {
		if (fabs(R(i,i)) > tol_) 
			is_ind[I(i)]=true;
	}
}

//compute the Cholesky decomposition of A = LL^T, L is lower triangular
void eigenCholDynamic( const int n, Real* A, Real* L);
void eigenCholFixed( const int n, Real* A, Real* L);
//call this one:
inline void chol( const int n, Real* A, Real* L) {
	eigenCholFixed(n,A,L);
}

#endif  //_WMRSIM_EIGENSOLVE_H_