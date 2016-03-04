#ifndef _WMRDE_OPTIONS_H_
#define _WMRDE_OPTIONS_H_

#include <string>

//OPTIONS, set these flags to 1 or 0
#define WMRDE_USE_QUATERNION 0 //else use Euler angles
#define WMRDE_USE_DOUBLE_PRECISION 1 //else use single (float) precision
#define WMRDE_ENABLE_ANIMATION 1 //include WmrAnimation, OGRE dependencies

//directory containing CAD model meshes of vehicles for animation
inline std::string meshDir() {
  return std::string("/home/neal/Projects/wmrde/CAD/");
}

//directory of other resources
inline std::string resourceDir()
{
  return std::string("/home/neal/Projects/wmrde/resource/");
}


//don't modify below here

//TODO, move this?
#if WMRDE_USE_QUATERNION
#define SIZEORIENT 4
#else
#define SIZEORIENT 3
#endif

//for output only in Debug mode
#ifdef _DEBUG
#define DEBUG_CERR(x) std::cerr << x << std::endl;
#else
#define DEBUG_CERR(x) //do {} while (0)
#endif


#endif  //_WMRDE_OPTIONS_H_
