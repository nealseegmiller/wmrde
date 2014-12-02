#include <demo/zoescene.h>


void zoeScene(const WmrModel& mdl, WmrAnimation& anim) {
	const bool draw_faces = false;
	const bool draw_edges = true;

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nw = mdl.get_nw();
	const int* wheelframeinds = mdl.get_wheelframeinds();

	const Frame* frames = mdl.get_frames();
	
	//add nodes for all frames
	for (int fi=0; fi < nf; fi++) {
		anim.addNode(frames[fi].parent_ind);
	}
	
	int i;
	std::string filename;

	i = mdl.nameToInd("Body");
	filename = VrmlDir() + std::string("Zoe/ZoeBody.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	i = mdl.nameToInd("SteerFront");
	filename = VrmlDir() + std::string("Zoe/ZoeFrontAxle.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	i = mdl.nameToInd("SteerRear");
	filename = VrmlDir() + std::string("Zoe/ZoeRearAxle.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	for (int wno = 0; wno < nw; wno++) {
		filename = VrmlDir() + std::string("Zoe/ZoeWheel.wrl");
		anim.addEntityVrml(anim.get_node(wheelframeinds[wno]), filename, -1, draw_faces, draw_edges);
	}

	//contact frames
	
	//add nodes for contact frames
	for (int wno = 0; wno < nw; wno++) {
		anim.addNode(wheelframeinds[wno]);
	}
	Real length = .25;

	Ogre::ColourValue C(0,1.0,0);

	//Ogre::Quaternion qrotx(Ogre::Degree( 90), Ogre::Vector3::UNIT_Y); //to make x axis arrow
	//Ogre::Quaternion qroty(Ogre::Degree(-90), Ogre::Vector3::UNIT_X); //to make y axis arrow

	for (int wno = 0; wno < nw; wno++) {
		//anim.addEntityArrow(anim.get_Node(nf+wno), length, .05*length, .25*length, .1*length, &qrotx, 0, &C);
		//anim.addEntityArrow(anim.get_Node(nf+wno), length, .05*length, .25*length, .1*length, &qroty, 0, &C);
		anim.addEntityArrow(anim.get_node(nf+wno), length, .05*length, .25*length, .1*length, 0, 0, &C);
	}

	//lines
	anim.addLine(1); //for body frame path
	anim.addLine(nw); //for contact point paths

}



