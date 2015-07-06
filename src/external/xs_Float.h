//http://stereopsis.com/sree/fpu2006.html
// ====================================================================================================================
// ====================================================================================================================
//  xs_Float.h
// ====================================================================================================================
// ====================================================================================================================
#ifndef _xs_FLOAT_H_
#define _xs_FLOAT_H_

#include "xs_Config.h"

// ====================================================================================================================
//  Defines
// ====================================================================================================================
#ifndef _xs_DEFAULT_CONVERSION
#define _xs_DEFAULT_CONVERSION      0
#endif //_xs_DEFAULT_CONVERSION


#if _xs_BigEndian_
	#define _xs_iexp_				0
	#define _xs_iman_				1
#else
	#define _xs_iexp_				1       //intel is little endian
	#define _xs_iman_				0
#endif //BigEndian_


#define _xs_doublecopysgn(a,b)      ((int32*)&a)[_xs_iexp_]&=~(((int32*)&b)[_xs_iexp_]&0x80000000) 
#define _xs_doubleisnegative(a)     ((((int32*)&a)[_xs_iexp_])|0x80000000) 

#ifndef _WIN32
#define __forceinline inline //to compile with GCC
#endif

// ====================================================================================================================
//  Constants
// ====================================================================================================================
const real64 _xs_doublemagic			= real64 (6755399441055744.0); 	    //2^52 * 1.5,  uses limited precisicion to floor
const real64 _xs_doublemagicdelta      	= (1.5e-8);                         //almost .5f = .5f + 1e^(number of exp bit)
const real64 _xs_doublemagicroundeps	= (.5f-_xs_doublemagicdelta);       //almost .5f = .5f - 1e^(number of exp bit)


// ====================================================================================================================
//  Prototypes
// ====================================================================================================================
int32 xs_CRoundToInt      (real64 val, real64 dmr =  _xs_doublemagic);
int32 xs_ToInt            (real64 val, real64 dme = -_xs_doublemagicroundeps);
int32 xs_FloorToInt       (real64 val, real64 dme =  _xs_doublemagicroundeps);
int32 xs_CeilToInt        (real64 val, real64 dme =  _xs_doublemagicroundeps);
int32 xs_RoundToInt       (real64 val);

//int32 versions
//commented out by Neal, caused linker errors
//int32 xs_CRoundToInt      (int32 val)   {return val;}
//int32 xs_ToInt            (int32 val)   {return val;}



// ====================================================================================================================
//  Fix Class
// ====================================================================================================================
template <int32 N> class xs_Fix
{
public:
    typedef int32 Fix;

    // ====================================================================================================================
    //  Basic Conversion from Numbers
    // ====================================================================================================================
    finline static Fix       ToFix       (int32 val)    {return val<<N;}
    finline static Fix       ToFix       (real64 val)   {return xs_ConvertToFixed(val);}

    // ====================================================================================================================
    //  Basic Conversion to Numbers
    // ====================================================================================================================
    finline static real64    ToReal      (Fix f)        {return real64(f)/real64(1<<N);}
    finline static int32     ToInt       (Fix f)        {return f>>N;}



protected:
    // ====================================================================================================================
    // Helper function - mainly to preserve _xs_DEFAULT_CONVERSION
    // ====================================================================================================================
    finline static int32 xs_ConvertToFixed (real64 val)
    {
    #if _xs_DEFAULT_CONVERSION==0
        return xs_CRoundToInt(val, _xs_doublemagic/(1<<N));
    #else
        return (long)((val)*(1<<N));
    #endif
    }
};





// ====================================================================================================================
// ====================================================================================================================
//  Inline implementation
// ====================================================================================================================
// ====================================================================================================================
finline int32 xs_CRoundToInt(real64 val, real64 dmr)
{
#if _xs_DEFAULT_CONVERSION==0
    val		= val + dmr;
	return ((int32*)&val)[_xs_iman_];
    //return 0;
#else
    return int32(floor(val+.5));
#endif
}


// ====================================================================================================================
finline int32 xs_ToInt(real64 val, real64 dme)
{
    /* unused - something else I tried...
            _xs_doublecopysgn(dme,val);
            return xs_CRoundToInt(val+dme);
            return 0;
    */

#if _xs_DEFAULT_CONVERSION==0
	return (val<0) ?   xs_CRoundToInt(val-dme) : 
					   xs_CRoundToInt(val+dme);
#else
    return int32(val);
#endif
}


// ====================================================================================================================
finline int32 xs_FloorToInt(real64 val, real64 dme)
{
#if _xs_DEFAULT_CONVERSION==0
    return xs_CRoundToInt (val - dme);
#else
    return floor(val);
#endif
}


// ====================================================================================================================
finline int32 xs_CeilToInt(real64 val, real64 dme)
{
#if _xs_DEFAULT_CONVERSION==0
    return xs_CRoundToInt (val + dme);
#else
    return ceil(val);
#endif
}


// ====================================================================================================================
finline int32 xs_RoundToInt(real64 val)
{
#if _xs_DEFAULT_CONVERSION==0
    return xs_CRoundToInt (val + _xs_doublemagicdelta);
#else
    return floor(val+.5);
#endif
}



// ====================================================================================================================
// ====================================================================================================================
#endif // _xs_FLOAT_H_
