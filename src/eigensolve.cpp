#include <wmrde/eigensolve.h>

void eigenSolveDynamic( const int nrows, const int ncols, Real* A, Real* b, Real* x ) {
	Eigen::Map<MatrixXr> A_(A,nrows,ncols);
	Eigen::Map<MatrixXr> b_(b,nrows,1);
	Eigen::Map<MatrixXr> x_(x,ncols,1);
	if ( nrows == ncols ) {
		x_ = A_.llt().solve(b_); //requires positive definite
	} else {
		x_ = A_.householderQr().solve(b_);
	}
}

void eigenSolveFixed( const int nrows, const int ncols, Real* A, Real* b, Real* x ) {

#ifdef FIXED_NROWS_0
	if (FIXED_NROWS_0 == nrows && FIXED_NCOLS_0 == ncols) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_0,FIXED_NCOLS_0>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_0,1>> b_(b);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NCOLS_0,1>> x_(x);
#if FIXED_NROWS_0 == FIXED_NCOLS_0
		x_ = A_.llt().solve(b_);
#else
		x_ = A_.householderQr().solve(b_);
#endif
		return;
	}
#endif


#ifdef FIXED_NROWS_1
	if (FIXED_NROWS_1 == nrows && FIXED_NCOLS_1 == ncols) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_1,FIXED_NCOLS_1>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_1,1>> b_(b);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NCOLS_1,1>> x_(x);
#if FIXED_NROWS_1 == FIXED_NCOLS_1
		x_ = A_.llt().solve(b_);
#else
		x_ = A_.householderQr().solve(b_);
#endif
		return;
	}
#endif


#ifdef FIXED_NROWS_2
	if (FIXED_NROWS_2 == nrows && FIXED_NCOLS_2 == ncols) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_2,FIXED_NCOLS_2>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_2,1>> b_(b);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NCOLS_2,1>> x_(x);
#if FIXED_NROWS_2 == FIXED_NCOLS_2
		x_ = A_.llt().solve(b_);
#else
		x_ = A_.householderQr().solve(b_);
#endif
		return;
	}
#endif


#ifdef FIXED_NROWS_3
	if (FIXED_NROWS_3 == nrows && FIXED_NCOLS_3 == ncols) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_3,FIXED_NCOLS_3>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_3,1>> b_(b);
		Eigen::Map<Eigen::Matrix<Real,FIXED_NCOLS_3,1>> x_(x);
#if FIXED_NROWS_3 == FIXED_NCOLS_3
		x_ = A_.llt().solve(b_);
#else
		x_ = A_.householderQr().solve(b_);
#endif
		return;
	}
#endif
	
#if PRINT_MATRIX_SIZE_IF_DYNAMIC
	std::cout << "resorting to eigenSolveDynamic, nrows = " << nrows << ", ncols = " << ncols << std::endl;
#endif

	eigenSolveDynamic(nrows, ncols, A, b, x); //fail safe
}


void eigenSubsetDynamic( const int nrows, const int ncols, Real* A, const Real tol, bool is_ind[]) {
	//obtain linearly independent subset of *rows* of A using QR decomposition

	Eigen::Map<MatrixXr> A_(A,nrows,ncols);
	Eigen::ColPivHouseholderQR<MatrixXr> QR(A_.transpose());

	MatrixXr R = QR.matrixQR().triangularView<Eigen::Upper>();
	Eigen::Matrix<int,Eigen::Dynamic,1> I = QR.colsPermutation().indices();

	eigenSetIsInd(R, I, tol, is_ind);
	
	//DEBUGGING
	//std::cout << "R=\n" << R << std::endl << std::endl;
	//std::cout << "I=\n" << I << std::endl;

	return;
}


//TODO

void eigenSubsetFixed( const int nrows, const int ncols, Real* A, const Real tol, bool is_ind[]) {
	//get linearly independent subset of *rows* of A using QR decomposition

#ifdef FIXED_NROWS_0
	if (FIXED_NROWS_0 == nrows && FIXED_NCOLS_0 == ncols) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_NROWS_0,FIXED_NCOLS_0>> A_(A);
		Eigen::ColPivHouseholderQR<Eigen::Matrix<Real,FIXED_NCOLS_0,FIXED_NROWS_0>> QR(A_.transpose());

		Eigen::Matrix<Real,FIXED_NCOLS_0,FIXED_NROWS_0> R = QR.matrixQR().triangularView<Eigen::Upper>();
		Eigen::Matrix<int,FIXED_NROWS_0,1> I = QR.colsPermutation().indices();
		
		eigenSetIsInd(R, I, tol, is_ind);

		return;
	}
#endif

#if PRINT_MATRIX_SIZE_IF_DYNAMIC
	std::cout << "resorting to eigenSubsetDynamic, nrows = " << nrows << ", ncols = " << ncols << std::endl;
#endif

	eigenSubsetDynamic( nrows, ncols, A, tol, is_ind ); //fail safe

	return;
}

//Cholesky decomposition
bool eigenCholDynamic( const int n, Real* A, Real* L) {
	Eigen::Map<MatrixXr> A_(A,n,n);
	Eigen::Map<MatrixXr> L_(L,n,n);
//	L_ = A_.llt().matrixL();

	Eigen::LLT<MatrixXr> A_llt(A_);
	L_ = A_llt.matrixL();

	return A_llt.info() != Eigen::Success;
}

bool eigenCholFixed( const int n, Real* A, Real* L) {
#ifdef FIXED_N_0
	if (FIXED_N_0 == n) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_N_0,FIXED_N_0>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_N_0,FIXED_N_0>> L_(L);
//		L_ = A_.llt().matrixL();

	  Eigen::LLT<MatrixXr> A_llt(A_);
	  L_ = A_llt.matrixL();

	  return A_llt.info() != Eigen::Success;
	}
#endif
#ifdef FIXED_N_1
	if (FIXED_N_1 == n) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_N_1,FIXED_N_1>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_N_1,FIXED_N_1>> L_(L);
//		L_ = A_.llt().matrixL();

	  Eigen::LLT<MatrixXr> A_llt(A_);
	  L_ = A_llt.matrixL();

	  return A_llt.info() != Eigen::Success;
	}
#endif
#ifdef FIXED_N_2
	if (FIXED_N_2 == n) {
		Eigen::Map<Eigen::Matrix<Real,FIXED_N_2,FIXED_N_2>> A_(A);
		Eigen::Map<Eigen::Matrix<Real,FIXED_N_2,FIXED_N_2>> L_(L);
//		L_ = A_.llt().matrixL();

	  Eigen::LLT<MatrixXr> A_llt(A_);
	  L_ = A_llt.matrixL();

	  return A_llt.info() != Eigen::Success;
	}
#endif

#if PRINT_MATRIX_SIZE_IF_DYNAMIC
	std::cout << "resorting to eigenCholDynamic, n = " << n << std::endl;
#endif

	return eigenCholDynamic(n,A,L); //fail safe

}
