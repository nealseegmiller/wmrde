//WmrAnimation.h
//This class is derived from the Ogre BaseApplication class used in Ogre tutorials
//This class is used to animate wheeled mobile robots simulations

#ifndef __WMRSIM_WMRANIMATION_H_
#define __WMRSIM_WMRANIMATION_H_

//OGRE
#include "animate/BaseApplication.h"
#include <OgreManualObject.h>
#include <OgreBillboardChain.h>
#include <OgreTechnique.h>
#include <animate/VrmlData.h>

//WmrSim
#include <contactgeom.h>
#include <surface/TriMeshSurf.h>
#include <surface/GridSurf.h>

class WmrAnimation : public BaseApplication
{
public:
	static const int MAXNN = 50; //max number of nodes
	static const int MAXNL = 5; //max number of lines
private:
	//TODO, move these members to local struct 'WmrScene'. one WmrScene object per WMR in scene
	int nn; //number of nodes (one per frame in WmrModel)
	Ogre::SceneNode* nodes[MAXNN];
	int parent_ind[MAXNN];

	//use BillBoardChain to draw thick lines. each line may have multiple chains
	Ogre::BillboardChain* lines[MAXNL];
	int nl; //number of lines

public:
    WmrAnimation(void);
    virtual ~WmrAnimation(void);
	virtual void go(void); //override go() in BaseApplication

	void start(void); //calls setup in BaseApplication, must call this first
	bool updateRender(void);

	//get
	Ogre::SceneNode* get_node(const int index) const { return nodes[index]; }
	const int get_nn() const { return nn; }
	const int get_nl() const { return nl; }
	
	//display WMR frames using VRML models
	void addNode(const int parent_node_index);
	void addEntityVrml(Ogre::SceneNode* node_to_attach_to, const std::string filename, const int flip_dim, const bool draw_faces, const bool draw_edges);

	//use arrows to display surface normal at contact point
	void addEntityBox(Ogre::SceneNode* node_to_attach_to, const Ogre::Vector3 L, 
		const Ogre::Quaternion* rot_off, const Ogre::Vector3* trans_off, const Ogre::ColourValue* C);
	void addEntityArrow(Ogre::SceneNode* NodeToAttachTo, const float length, const float width, const float hlength, const float hwidth, 
		const Ogre::Quaternion* rot_off, const Ogre::Vector3* trans_off, const Ogre::ColourValue* C);
	//display terrain surface
	void addEntitySurface( Surface* surf );
	void addEntityTriMeshSurf( TriMeshSurf* surf );
	void addEntityGridSurf( GridSurf* surf );

	void updateNodeTransform(const int node_index, const HomogeneousTransform HT) {
		Ogre::Quaternion quat( Ogre::Matrix3( 
			HT[0], HT[0+COL1], HT[0+COL2], 
			HT[1], HT[1+COL1], HT[1+COL2],
			HT[2], HT[2+COL1], HT[2+COL2] ) );
		Ogre::Vector3 pos( HT[0+COL3], HT[1+COL3], HT[2+COL3] );

		nodes[node_index]->setOrientation( quat );
		nodes[node_index]->setPosition( pos );
	}

	//use lines to display the trajectory of frames
	void addLine( const int num_chains );
	void addChainElement( const int lineno, const int chainno, const Ogre::Vector3 pt, const float thickness, Ogre::ColourValue* C);

	//updates nodes and lines
	void updateNodesLines( const int nf, const HomogeneousTransform HT_to_parent[], const int num_contacts, const ContactGeom* contacts );

protected:
    virtual void createScene(void);
};

#endif
