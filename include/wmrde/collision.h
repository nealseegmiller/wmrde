//collision.h
//functions for collision detection between wheels/tracks and terrain

#ifndef _WMRDE_COLLISION_H_
#define _WMRDE_COLLISION_H_

#include <assert.h>

#include <wmrde/surface/Surface.h>
#include <wmrde/contactgeom.h>
#include <wmrde/linesearch.h>

void updateWheelContactGeomDiscretize(const SurfaceVector& surfaces, const HomogeneousTransform HT_wheel_to_world, const Real radius,  //input
	WheelContactGeom& contact); //output

void updateWheelContactGeomRoot(const SurfaceVector& surfaces, const HomogeneousTransform HT_wheel_to_world, const Real radius, 
	WheelContactGeom& contact);

inline void updateWheelContactGeom(const SurfaceVector& surfaces, const HomogeneousTransform HT_wheel_to_world, const Real radius, 
	WheelContactGeom& contact) {
	//uncomment one of the following:
//	updateWheelContactGeomDiscretize(surfaces, HT_wheel_to_world, radius, contact); //matches MATLAB
	updateWheelContactGeomRoot(surfaces, HT_wheel_to_world, radius, contact); //faster!
}

//convert contact angle to point in wheel coordinates
//convention:
//contact angle=0 corresponds to cp=[0,0,-rad]'
//ccw about y axis is positive
inline void contactAngleToPoint(const Real radius, const Real angle, Vec3 point_wheel) {
	setVec3(-sin(angle)*radius, 0, -cos(angle)*radius, point_wheel);
}

void HTContactToWheel(const Vec3 pt, const Vec3 N, HomogeneousTransform HT_contact_to_wheel);

void whichPointsInContactWheel(const int min_npic, const int nw, WheelContactGeom* contacts);
void whichPointsInContactTrack(const int min_npic, const int nt, TrackContactGeom* contacts);
void whichPointsInContact(const int min_npic, const int np, const Vec3 pts_world[], const Real dz[], int incontact[]);

//for tracks
void initTrackContactGeom(const Real rad, const Real rad2, const Real L, TrackContactGeom& contact);
int HTContactToTrack(const Real rad, const Real rad2, const Real L, const int npflat, const bool sides[], HomogeneousTransform HT_track[]);
void updateTrackContactGeom(const SurfaceVector& surfaces, const HomogeneousTransform HT_track_to_world, //input
	TrackContactGeom& contact); //output

#endif //_WMRDE_COLLISION_H_
