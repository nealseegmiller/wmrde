//index.h
//implements built-in MATLAB functionality for indexing arrays

#ifndef _WMRSIM_INDEX_H_
#define _WMRSIM_INDEX_H_

#include <algorithm>
#include <common/common.h>

//inline functions

//val:	n evenly distributed values from lo to hi
template<typename Type>
inline void linspace(const Type lo, const Type hi, const int n, Type val[]) {
	Type d = (hi-lo)/(Type(n)-1);
	val[0]=lo;
	for (int i=1; i<n; i++)
		val[i] = val[i-1]+d; 
}

//(return) index of minimum value in array
template<typename Type>
inline int findMin(const int n, const Type val[]) {
	Type min_val=val[0];
	int min_idx=0;
	for (int i=1; i<n; i++) {
		if (val[i] < min_val) {
			min_val=val[i];
			min_idx=i;
		}
	}
	return min_idx;
}

template<typename Type>
inline int findMax(const int n, const Type val[]) {
	Type max_val=val[0];
	int max_idx=0;
	for (int i=1; i<n; i++) {
		if (val[i] > max_val) {
			max_val=val[i];
			max_idx=i;
		}
	}
	return max_idx;
}

//(return)	index of minimum *absolute* value in array
template<typename Type>
inline int findMinAbs(const int n, const Type val[]) {
	Type val_;
	Type min_val=val[0];
	int min_idx=0;
	for (int i=1; i<n; i++) {
		val_ = TYPESIGN(val[i],Type)*val[i]; //absolute value
		if (val_ < min_val) {
			min_val=val_;
			min_idx=i;
		}
	}
	return min_idx;
}

template<typename Type>
inline int findMaxAbs(const int n, const Type val[]) {
	Type val_;
	Type max_val=val[0];
	int max_idx=0;
	for (int i=1; i<n; i++) {
		val_ = TYPESIGN(val[i],Type)*val[i]; //absolute value
		if (val_ > max_val) {
			max_val=val_;
			max_idx=i;
		}
	}
	return max_idx;
}

//logical indexing
//(return) number of true elements in logical array
inline int logicalCount( const int n, const bool logical[] ) {
	int ntrue=0;
	for (int i=0; i<n; i++) {
		ntrue += logical[i];
	}
	return ntrue;
}

//inds:		indices of true elements in logical array
//(return)	number of true elements
inline int logicalFind( const int n, const bool logical[], int inds[] ) {
	int ntrue=0;
	for (int i=0; i<n; i++) {
		if (logical[i]) {
			inds[ntrue]=i;
			ntrue++;
		}
	}
	return ntrue;
}

//use logical to index array_in, copy to array_out
template<typename Type>
inline int logicalIndexIn( const int n, const bool logical[], const Type array_in[], Type array_out[] ) {
	int ntrue=0;
	for (int i=0; i<n; i++) {
		if (logical[i]) {
			array_out[ntrue]=array_in[i];
			ntrue++;
		}
	}
	return ntrue;
}

//use logical to index array_out, copy from array in
template<typename Type>
inline int logicalIndexOut( const int n, const bool logical[], const Type array_in[], Type array_out[] ) {
	int ntrue=0;
	for (int i=0; i<n; i++) {
		if (logical[i]) {
			array_out[i]=array_in[ntrue];
			ntrue++;
		}
	}
	return ntrue;
}

template <typename Type>
inline void sortIndex(const int n, const Type val[], int idx[], Type val_sort[]) {
	//initialize
	for (int i=0; i<n; i++)
		idx[i] = i;

	std::sort(idx, idx + n, [&val](int i1, int i2) { return val[i1] < val[i2]; } );

	if (val_sort != 0) { //not null
		for (int i=0; i<n; i++)
			val_sort[i] = val[idx[i]];
	}
}


#endif