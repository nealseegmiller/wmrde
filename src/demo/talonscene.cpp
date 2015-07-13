#include <wmrde/demo/talonscene.h>


void talonScene(const WmrModel& mdl, const TrackContactGeom contacts[], WmrAnimation& anim) {
	const bool draw_faces = false;
	const bool draw_edges = true;

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nt = mdl.get_nt();
	const int* sprocketframeinds = mdl.get_sprocketframeinds();

	const Frame* frames = mdl.get_frames();
	
	//add nodes for all frames
	for (int fi=0; fi < nf; fi++) {
		anim.addNode(frames[fi].parent_ind);
	}
	
	int i;
	std::string filename;

	i = mdl.nameToInd("Body");
	filename = CADdir() + std::string("Talon/TalonBody.wrl");
	anim.addEntityVrml(anim.get_node(i), filename, -1, draw_faces, draw_edges);

	for (int tno = 0; tno < nt; tno++) {
		filename = CADdir() + std::string("Talon/TalonSprocket.wrl");
		anim.addEntityVrml(anim.get_node(sprocketframeinds[tno]), filename, -1, draw_faces, draw_edges);
	}

	//for arrows
	Ogre::ColourValue C(0,1.0,0);
	Real length = .1;

	int node_idx = nf-1; //node index

	for (int tno = 0; tno < nt; tno++) {
		int fi = sprocketframeinds[tno];

		//add node for track frame
		anim.addNode(frames[fi].parent_ind);
		node_idx++;

		int track_node_idx = node_idx;

		anim.updateNodeTransform(node_idx,frames[fi].HT_parent_jd0);

		filename = CADdir() + std::string("Talon/TalonTrack.wrl");
		anim.addEntityVrml(anim.get_node(node_idx), filename, -1, draw_faces, draw_edges);

		int np = contacts[tno].get_np();
		for (int pno = 0; pno < np; pno++) {

			//add node for contact frame
			anim.addNode(track_node_idx);
			node_idx++;

			anim.updateNodeTransform(node_idx,contacts[tno].HT_track[pno]);

			//z axis arrow
			anim.addEntityArrow(anim.get_node(node_idx), length, .05*length, .25*length, .1*length, 0, 0, &C);
		}

	}
	
	//lines
	anim.addLine(1); //for body frame path

}
