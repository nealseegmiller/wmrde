#include <wmrde/algebra/matrix.h>

//solve a linear system given the Cholesky decomposition
//A*x=b
//A=L*L'	L is lower triangular
//L*y=b		solve for y by forward substitution
//L'*x=y	solve for x by back substitution
void cholSolve(const int n, const Real L[], const Real b[], Real x[]) {

	//forward substitution, L*y=b
	//store y in x
	for (int i=0; i<n; i++) {
		x[i] = b[i];
		//TODO, not contiguous, can't vectorize
		for (int j=0; j<i; j++) {
			x[i] -= L[S2I(i,j,n)]*x[j];
		}
		x[i] /= L[S2I(i,i,n)];
	}

	//std::cout << "b=\n"; printMatReal(n,1,b,-1,-1);
	//std::cout << "y=\n"; printMatReal(n,1,x,-1,-1);

	//back substitution, L'*x=y
	for (int i=n-1; i>=0; i--) {
		//x[i] = x[i];
		for (int j=i+1; j<n; j++) {
			x[i] -= L[S2I(j,i,n)]*x[j]; //transpose
		}
		x[i] /= L[S2I(i,i,n)];
	}

	//std::cout << "x=\n"; printMatReal(n,1,x,-1,-1);

	return;
}

//solve linear system given the Cholesky decomposition, where solution is a matrix
//A*X=B, B is n x ncols
//A=L*L'
//y is an intermediate variable
void cholSolveMat(const int n, const Real L[], const int ncols, const Real B[], Real X[]) {
	for (int i=0; i<ncols; i++) 
		cholSolve(n,L,B+S2I(0,i,n),X+S2I(0,i,n));
}

//print a matrix of real numbers
//nrows:	number of rows
//ncols:	number of cols
//M:		nrows*ncols array of values
//precision:	number of digits after decimal, -1 to use default
//width:		-1 to use default, 0 to automatically set min width
void printMatReal(const int nrows, const int ncols, const Real M[], int precision, int width) {

	//set to default
	if (precision < 0) { precision = 4; }
	if (width < 0) { width = 0; }
	
	//backup format
	std::ios_base::fmtflags old_flags = std::cout.flags();
	std::streamsize old_precision = std::cout.precision();

	//set format
	std::cout << std::fixed << std::setprecision(precision);

	if (width == 0) {
		Real max;
		max = M[findMaxAbs(nrows*ncols,M)];
		max = fabs(max);
		int lmax = (int) log10(max); //truncates toward zero
		if (max > 0 && (lmax < -precision || lmax > 10)) {
			//numbers are nonzero and very small or large, so use scientific notation
			std::cout << std::scientific;
			//negative sign + single digit before decimal + decimal point + digits after decimal + exponent
			width = 3 + precision + 5;
		} else {
			//negative sign + digits before decimal + digits after decimal
			width = 1 + (std::max(lmax,0)+1) + precision;
			if (precision > 0) 
				width++; //+ decimal point
		}
	}

	//print
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++) {
			if ( fabs( M[i + j*nrows] ) > 0.0 ) {
				std::cout << std::setw(width) << M[S2I(i,j,nrows)] << ' ';
			} else {
				std::cout << std::setw(width) << (int) 0 << ' ';
			}
		}
		std::cout << std::endl;
	}
	//reset format
	std::cout.flags(old_flags);
	std::cout.precision(old_precision);
}

//print a matrix of integers
void printMatInt(const int nrows, const int ncols, const int M[], int width) {

	//set to default
	if (width < 0) { width = 0; }

	if (width == 0) {
		int max;
		max = M[findMaxAbs(nrows*ncols,M)];
		max = abs(max);
		int lmax = (int) log10((Real) max); //truncates toward zero
		width = 1 + (std::max(lmax,0)+1); //negative sign + number of digits
	}

	//print
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++) {
			std::cout << std::setw(width) << M[S2I(i,j,nrows)] << ' ';
		}
		std::cout << std::endl;
	}
}

//print a matrix of boolean values
void printMatBool(const int nrows, const int ncols, const bool M[], int width) {

	//set to default
	if (width < 1) width = 1;
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++) {
			std::cout << std::setw(width) << (int) M[S2I(i,j,nrows)] << ' ';
		}
		std::cout << std::endl;
	}

}
