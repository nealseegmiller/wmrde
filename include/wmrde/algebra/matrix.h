//matrix.h
//linear algebra for vectors/matrices of any size
//these are faster than using Eigen dynamic matrices for small matrices. Tested 10x10, 20x20
//auto-vectorized where possible. Visual Studio C++ command line /Qvec-report:2
//Visual Studio vectorizer messages: 
//http://msdn.microsoft.com/en-us/library/jj658585.aspx
#ifndef WMRDE_MATRIX_H_
#define WMRDE_MATRIX_H_

//for printMat*
#include <iostream>
#include <iomanip>
#include <cmath>

//#include <wmrde/options.h>
#include <wmrde/common.h>
#include <wmrde/util/index_util.h> //for findMaxAbs() for print functions

namespace wmrde
{

//column-major order assumed!
//subscripts to index
#define S2I(ri,ci,nrows) (ri+(ci)*(nrows))

//vectors are 1D arrays, vector operations are equivalent to column operations with ci=0


//SET FUNCTIONS
//all set operations, loop not vectorized due to reason 1300 (no computation)
template<typename Type>
inline void setVec(const int n, const Type val, Type Dest[]) {
	for (int i=0; i < n; i++) 
		Dest[i] = val;
}

template<typename Type>
inline void setMat(const int nrows, const int ncols, const Type val, Type Dest[] ) { 
	int n = nrows*ncols;
	for (int i=0; i < n; i++) 
		Dest[i] = val;
}

//set row is faster than set column, why?
template<typename Type>
inline void setMatRow(const int nrows, const int ncols, const int ri, const Type val, Type Dest[] ) { 
	for (int j=0; j < ncols; j++) 
		Dest[S2I(ri,j,nrows)] = val;
}

template<typename Type>
inline void setMatCol(const int nrows, const int ci, const Type val, Type Dest[] ) { 
	for (int i=0; i < nrows; i++) 
		Dest[S2I(i,ci,nrows)] = val;
}

template<typename Type>
inline void setMatBlock(const int nrows, const int ri, const int ci, const int nrows_block, const int ncols_block, const Type val, Type Dest[] ) { 
	for (int i=0; i < nrows_block; i++) {
		for (int j=0; j < ncols_block; j++) {
			Dest[S2I(i+ri,j+ci,nrows)] = val;
		}
	}
}



//COPY FUNCTIONS
//avoid copying when possible, faster to write directly to destination
template<typename Type>
inline void copyVec(const int n, const Type Source[], Type Dest[]) { 
	for (int i=0; i < n; i++) 
		Dest[i] = Source[i];
}

template<typename Type>
inline void copyMat(const int nrows, const int ncols, const Type Source[], Type Dest[]) {

	//fastest
	//int n=nrows*ncols; //using intermediate var is much slower, why?
	for (int i=0; i < nrows*ncols; i++) 
		Dest[i] = Source[i];

}

template<typename Type>
inline void copyMatRow( const int nrows_source, const int ncols, const int ri_source, const Type Source[], const int nrows_dest, const int ri_dest, Type Dest[] ) {
	for (int j=0; j < ncols; j++) 
		Dest[S2I(ri_dest,j,nrows_dest)] = Source[S2I(ri_source,j,nrows_source)];
}

//copy col is faster than copy row, why?
template<typename Type>
inline void copyMatCol(const int nrows, const int ci_source, const Type Source[], const int ci_dest, Type Dest[]) { 
	for (int i=0; i < nrows; i++) 
		Dest[S2I(i,ci_dest,nrows)] = Source[S2I(i,ci_source,nrows)];

}

template<typename Type>
inline void copyMatBlock(const int nrows_source, const int ri_source, const int ci_source, const int nrows_block, const int ncols_block, const Type Source[], 
				  const int nrows_dest, const int ri_dest, const int ci_dest, Type Dest[]) {

	//faster, why?
	for (int j=0; j < ncols_block; j++) { //col index
		for (int i=0; i < nrows_block; i++) { //row index
			Dest[S2I(i+ri_dest,j+ci_dest,nrows_dest)] = Source[S2I(i+ri_source,j+ci_source,nrows_source)];
		}
	}

	//for (int j=0; j < ncols_block; j++) { //col index
	//	copyVec( nrows_block, Source+S2I(ri_source,j+ci_source,nrows_source), Dest+S2I(ri_dest,j+ci_dest,nrows_dest) );
	//}
}


//ADD FUNCTIONS
//same parameters as copy functions, plus a coefficient (m = -1 for subtraction)
template<typename Type>
inline void addmVec(const int n, const Type Source[], const Type m, Type Dest[]) { 
	for (int i=0; i < n; i++) 
		Dest[i] += m*Source[i];
}

template<typename Type>
inline void addmMat(const int nrows, const int ncols, const Type Source[], const Type m, Type Dest[]) {
	//fastest
	int n = nrows*ncols; //much slower without this
	for (int i=0; i < n; i++) 
		Dest[i] += m*Source[i]; //vectorized

}

template<typename Type>
inline void addmMatRow( const int nrows_source, const int ncols, const int ri_source, const Type Source[], const int nrows_dest, const int ri_dest, const Type m, Type Dest[] ) {
	for (int j=0; j < ncols; j++) 
		Dest[S2I(ri_dest,j,nrows_dest)] += m*Source[S2I(ri_source,j,nrows_source)];
}

template<typename Type>
inline void addmMatCol(const int nrows, const int ci_source, const Type Source[], const int ci_dest, const Type m, Type Dest[]) { 
	//for (int i=0; i < nrows; i++) 
	//	Dest[S2I(i,ci_dest,nrows)] += m*Source[S2I(i,ci_source,nrows)];

	//vectorized, faster
	const Type* source_ptr = Source+S2I(0,ci_source,nrows);
	Type* dest_ptr = Dest+S2I(0,ci_dest,nrows);
	addmVec(nrows, source_ptr, m, dest_ptr);
	
}

template<typename Type>
inline void addmMatBlock(const int nrows_source, const int ri_source, const int ci_source, const int nrows_block, const int ncols_block, const Type Source[], 
				  const int nrows_dest, const int ri_dest, const int ci_dest, const Type m, Type Dest[]) {

	//to vectorize, if sufficient loop iterations
	const Type* source_ptr;
	Type* dest_ptr;
	for (int j=0; j < ncols_block; j++) {
		source_ptr = Source + S2I(ri_source,j+ci_source,nrows_source);
		dest_ptr = Dest + S2I(ri_dest,j+ci_dest,nrows_dest);
		addmVec( nrows_block, source_ptr, m, dest_ptr );
	}
}



//MULTIPLY BY CONSTANT FUNCTIONS
//same parameters as set functions
template<typename Type>
inline void mulcVec(const int n, const Type m, Type Dest[]) {
	for (int i=0; i < n; i++) 
		Dest[i] *= m;
}

template<typename Type>
inline void mulcMat(const int nrows, const int ncols, const Type m, Type Dest[] ) { 
	
	//int n = nrows*ncols; //much slower without this
	//for (int i=0; i < n; i++) 
	//	Dest[i] *= m; //vectorized, why slow?

	//fastest, why?
	for (int j=0; j < ncols; j++) {
		for (int i=0; i < nrows; i++) {
			Dest[S2I(i,j,nrows)] *= m;
		}
	}

}

template<typename Type>
inline void mulcMatRow(const int nrows, const int ncols, const int ri, const Type m, Type Dest[] ) { 
	for (int j=0; j < ncols; j++) 
		Dest[S2I(ri,j,nrows)] *= m;
}

template<typename Type>
inline void mulcMatCol(const int nrows, const int ci, const Type m, Type Dest[] ) { 
	for (int i=0; i < nrows; i++) 
		Dest[S2I(i,ci,nrows)] *= m; //not vectorized

	//mulcVec( nrows, m, Dest+S2I(0,ci,nrows) ); //vectorized, why slower?
}

template<typename Type>
inline void mulcMatBlock(const int nrows, const int ri, const int ci, const int nrows_block, const int ncols_block, const Type m, Type Dest[] ) { 

	//fastest, why?
	for (int j=0; j < ncols_block; j++) {
		for (int i=0; i < nrows_block; i++) {
			Dest[S2I(i+ri,j+ci,nrows)] *= m;
		}
	}

	//for (int j=0; j < ncols_block; j++) {
	//	mulcVec( nrows_block, m, Dest+S2I(ri,j+ci,nrows) ); //to vectorize
	//}
}




//COPY TRANSPOSE OPERATIONS

template<typename Type>
inline void copyTMat(const int nrows, const int ncols, const Type Source[], Type Dest[]) {
	//TODO, fastest way?
	for (int i=0; i < nrows; i++) {
		for (int j=0; j < ncols; j++) {
			Dest[S2I(j,i,ncols)] = Source[S2I(i,j,nrows)];
		}
	}
}

template<typename Type>
inline void copyRowToVec( const int nrows, const int ncols, const int ri, const Type Source[], Type Dest[] ) {
	for (int i=0; i < ncols; i++) 
		Dest[i] = Source[S2I(ri,i,nrows)];
}

template<typename Type>
inline void copyVecToRow( const Type Source[], const int nrows, const int ncols, const int ri, Type Dest[] ) {
	for (int i=0; i < ncols; i++) 
		Dest[S2I(ri,i,nrows)] = Source[i];
}





//DOT PRODUCT OPERATIONS
//a and b are both vectors
template<typename Type>
inline Type dotVec(const int n, const Type a[], const Type b[]) {
	Type p=0; //product
	for (int i=0; i<n; i++) 
		p += a[i]*b[i];
	return p;
}


//MATRIX MULTIPLICATION OPERATIONS
//TODO, remove coefficient m?

//Ab = A*b
template<typename Type>
inline void multMatVec(const int nrows, const int ncols, const Type A[], const Type b[], const Type m, Type Ab[]) {
	//for (int i=0; i < nrows; i++) { //row in A
	//	Ab[i] = m*dotRowVec( ncols, nrows, A+S2I(i,0,nrows), b );
	//}

	//to vectorize
	setVec(nrows, 0.0, Ab);
	const Type* source_ptr;
	for (int j=0; j < ncols; j++) { //col in A
		source_ptr = A+S2I(0,j,nrows);
		addmVec(nrows, source_ptr, m*b[j], Ab);
	}
}

//ATb = (A^T)*b
template<typename Type>
inline void multMatTVec(const int nrows, const int ncols, const Type A[], const Type b[], const Type m, Type ATb[]) {
	const Type* ptr; //necessary to vectorize
	for (int j=0; j < ncols; j++) { //col in A
		ptr = A+S2I(0,j,nrows);
		ATb[j] = m*dotVec( nrows, ptr, b ); //vectorized

		//ATb[j] = m*dotVec( nrows, A+S2I(0,j,nrows), b );
	}
}

//AB = A*B
//A is nrows_a x ncols_a
//B is ncols_a x ncols_b
//AB is nrows_a x ncols_b
template<typename Type>
inline void multMatMat(const int nrows_a, const int ncols_a, const Type A[], const int ncols_b, const Type B[], const Type m, Type AB[]) {
	for (int j=0; j < ncols_b; j++) //col in B, AB
		multMatVec( nrows_a, ncols_a, A, B+S2I(0,j,ncols_a), m, AB+S2I(0,j,nrows_a) );
	
	//pointers not necessary to vectorize?

}


//ATB = (A^T)*B
//A is nrows_a x ncols_a
//B is nrows_a x ncols_b
//ATB is ncols_a x ncols_b
template<typename Type>
inline void multMatTMat(const int nrows_a, const int ncols_a, const Type A[], const int ncols_b, const Type B[], const Type m, Type ATB[]) {
	const Type* B_ptr; //necessary to vectorize
	for (int j=0; j < ncols_b; j++) { //col in B, ATB
		B_ptr = B+S2I(0,j,nrows_a);
		multMatTVec( nrows_a, ncols_a, A, B_ptr, m, ATB+S2I(0,j,ncols_a) ); //vectorized

		//multMatTVec( nrows_a, ncols_a, A, B+S2I(0,j,nrows_a), m, ATB+S2I(0,j,ncols_a) );
	}
}

void cholSolve(const int n, const Real L[], const Real b[], Real x[]);
void cholSolveMat(const int n, const Real L[], const int ncols, const Real B[], Real X[]);

void printMatReal(const int nrows, const int ncols, const Real M[], int precision, int width );
void printMatInt(const int nrows, const int ncols, const int M[], int width);
void printMatBool(const int nrows, const int ncols, const bool M[], int width);

} //namespace

#endif  //WMRDE_MATRIX_H_
