#include <wmrde/demo/rockyscene.h>


void rockyScene(const WmrModel& mdl, WmrAnimation& anim) {
	bool draw_faces = false;
	bool draw_edges = true;

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
	filename =  CADdir() + std::string("Rocky7/Rocky7Body.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	//rockers
	i = mdl.nameToInd("D1");
	filename = CADdir() + std::string("Rocky7/Rocky7Rocker.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	i = mdl.nameToInd("D2");
	filename = CADdir() + std::string("Rocky7/Rocky7Rocker.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, 1, draw_faces, draw_edges);

	//bogies
	i = mdl.nameToInd("B1");
	filename = CADdir() + std::string("Rocky7/Rocky7Bogie.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	i = mdl.nameToInd("B2");
	filename = CADdir() + std::string("Rocky7/Rocky7Bogie.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, 1, draw_faces, draw_edges);

	//steering brackets
	i = mdl.nameToInd("S1");
	filename = CADdir() + std::string("Rocky7/Rocky7Bracket.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	i = mdl.nameToInd("S2");
	filename = CADdir() + std::string("Rocky7/Rocky7Bracket.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, 1, draw_faces, draw_edges);

	//wheels
	filename = CADdir() + std::string("Rocky7/Rocky7Wheel.wrl");
	//left wheels: A1, A3, A5
	for (int wno=0; wno<nw; wno=wno+2) {
		i = wheelframeinds[wno];
		anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);
	}

	//right wheels: A2, A4, A6
	for (int wno=1; wno<nw; wno=wno+2) {
		i = wheelframeinds[wno];
		anim.addEntityVrml(anim.get_node(i), filename, 1, draw_faces, draw_edges);
	}

	//contact frames
	
	//add nodes for contact frames
	for (int wno = 0; wno < nw; wno++) {
		anim.addNode(wheelframeinds[wno]);
	}
	Real length = .1;
	Ogre::ColourValue C(0,1.0,0);
	for (int wno = 0; wno < nw; wno++) {
		anim.addEntityArrow(anim.get_node(nf+wno), length, .05*length, .25*length, .1*length, 0, 0, &C);
	}

	//lines
	anim.addLine(1); //for body frame path
	anim.addLine(nw); //for contact point paths

}

