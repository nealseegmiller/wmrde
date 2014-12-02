//contactgeom.h
//WheelContactGeom and TrackContactGeom classes
//required inputs to kinematics and dynamics functions for contact constraints

#ifndef _WMRSIM_CONTACTGEOM_H_
#define _WMRSIM_CONTACTGEOM_H_

#include <algebra/transform.h>

//base class
//avoid conflict with Open Dynamics Engine dContactGeom
class ContactGeom {
public:
	static const int MAXNP = 11; //max number of contact points

	virtual int get_np() const = 0; //get number of contact points
};


class WheelContactGeom : public ContactGeom {
public:
	Real dz; //contact height error
	Real angle; //contact angle
	HomogeneousTransform HT_wheel; //transform from contact to wheel coords
	HomogeneousTransform HT_world; //transform from contact to world coords
	bool incontact;

	//methods
	int get_np() const { return 1; }
};

class TrackContactGeom : public ContactGeom {
private:
	int np;		//number of contact points

public:
	Real dz[MAXNP];
	HomogeneousTransform HT_track[MAXNP];	//transform from contact to track coords
	HomogeneousTransform HT_world[MAXNP];	//transform from contact to world coords
	bool incontact[MAXNP];

	//methods
	void set_np(int number_of_points) { np = number_of_points; }

	int get_np() const { return np; }

	
};

#endif