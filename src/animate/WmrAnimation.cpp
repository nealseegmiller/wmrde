#include "animate/WmrAnimation.h"

//-------------------------------------------------------------------------------------
WmrAnimation::WmrAnimation(void)
{
	nn=0;
}
//-------------------------------------------------------------------------------------
WmrAnimation::~WmrAnimation(void)
{
}

//-------------------------------------------------------------------------------------
void WmrAnimation::createScene(void)
{

	// Set the scene's ambient light
    mSceneMgr->setAmbientLight(Ogre::ColourValue(0.5f, 0.5f, 0.5f));

	// Create a Light and set its position
    Ogre::Light* light = mSceneMgr->createLight("MainLight");
    light->setPosition(20.0f, 80.0f, 50.0f);

}

void WmrAnimation::addNode(const int parent_node_index) {
	assert(nn < MAXNN);

	if (parent_node_index > -1) {
		nodes[nn] = nodes[parent_node_index]->createChildSceneNode();
	} else { //body frame
		nodes[nn] = mSceneMgr->getRootSceneNode()->createChildSceneNode();
	}
	parent_ind[nn] = parent_node_index;
	nn++;
}

//read in VRML (.wrl) file and attach to node
//node_to_attach_to:	a pointer to a previously added node
//filename:				.wrl file name
//flip_dim:				0-2, dimension to flip coordinates about. no flip if < 0
//draw_faces:			
//draw_edges:
//TODO, rot_off, trans_off inputs
void WmrAnimation::addEntityVrml(Ogre::SceneNode* node_to_attach_to, const std::string filename, const int flip_dim, const bool draw_faces, const bool draw_edges) {

	//load vrml
	VrmlData<float> vrml;
	int success = vrml.readfile(filename.c_str());

	Ogre::ManualObject* manual = mSceneMgr->createManualObject();
		
	//duplicate vertices for multiple normal directions, better lighting
	//and for color per triangle vs. per vertex

	int np = vrml.nparts();
	for (int part_no=0; part_no<np; part_no++) { //loop over parts
			
		manual->setDynamic(true);

		int nt = vrml.ntriangles(part_no); //number of triangles

		if (draw_faces) {
			//draw faces

			manual->estimateVertexCount(nt*3);
			//manual->estimateIndexCount(nt*3);
			manual->estimateIndexCount(0);

			manual->begin("MyMaterials/FlatVertexColor", Ogre::RenderOperation::OT_TRIANGLE_LIST);

			for (int tno=0; tno < nt; tno++) { //loop over triangles
				
				float normal[3];
				float color[3]; //rgb

				vrml.getTriangleNormal(part_no,tno,normal[0],normal[1],normal[2]);
				vrml.getTriangleColor(part_no,tno,color[0],color[1],color[2]);

				if (flip_dim > 0) {
					//flip normal
					normal[0] = -normal[0];
					normal[1] = -normal[1];
					normal[2] = -normal[2];
				}
				
				
				for (int vno=0; vno < 3; vno++) { //loop over vertices
					float pos[3];

					vrml.getTriangleVertex(part_no,tno,vno,pos[0],pos[1],pos[2]); 

					if (flip_dim > 0)
						pos[flip_dim] = -pos[flip_dim];

					manual->position(pos[0],pos[1],pos[2]);
					manual->normal(normal[0],normal[1],normal[2]);
					manual->colour(color[0],color[1],color[2]);
				}
			}

			//not necessary?
			for (int tno=0; tno < nt; tno++) {
				int vi = tno*3; //vertex index in the entire list
				if (flip_dim > 0) {
					//swap indices for anti-clockwise order
					manual->triangle(vi, vi+2, vi+1); 
				} else {
					manual->triangle(vi, vi+1, vi+2); 
				}
			}

			manual->end();
		}

		//TODO, remove this?
		if (draw_edges) {
			//draw edges
			//doesn't look good with faces
			float tol = Ogre::Math::DegreesToRadians(15);
			int nv = vrml.ncoords(part_no);
			//normal indices for triangles that contain edge
			std::vector<int> normal_idx_1(nv*nv,-1); //1st
			std::vector<int> normal_idx_2(nv*nv,-1); //2nd
			//vertex indices for edges to draw
			std::vector<int> evert_idx_1; //1st vertex index
			std::vector<int> evert_idx_2; //2nd
				
			for (int tno=0; tno < nt; tno++) { //loop over triangles
				//loop over triangle edges
				for (int vno_1 = 0; vno_1 < 3; vno_1++) { //1st vertex no (in the triangle)
					int vno_2; //2nd vertex no
					if (vno_1 < 2) 
						vno_2 = vno_1 + 1; 
					else
						vno_2 = 0;

					//vertex indices in the entire list
					int v1 = vrml.getTriangleVertexIndex(part_no,tno,vno_1);
					int v2 = vrml.getTriangleVertexIndex(part_no,tno,vno_2);

					//sort
					if (v1>v2) 
						std::swap(v1,v2);

					if ( normal_idx_1[v1 + v2*nv] == -1 ) {
						normal_idx_1[v1 + v2*nv] = vrml.getTriangleNormalIndex(part_no, tno);
					} else {
						normal_idx_2[v1 + v2*nv] = vrml.getTriangleNormalIndex(part_no, tno);

						float N1[3], N2[3]; //normal vectors
						vrml.getNormal(part_no, normal_idx_1[v1 + v2*nv],N1[0],N1[1],N1[2]);
						vrml.getNormal(part_no, normal_idx_2[v1 + v2*nv],N2[0],N2[1],N2[2]);

						float angle = acos(N1[0]*N2[0] + N1[1]*N2[1] + N1[2]*N2[2]);
						if (fabs(angle) > tol) {
							evert_idx_1.push_back(v1);
							evert_idx_2.push_back(v2);
						}
					}
				}
			}

			Ogre::ColourValue C(1.0,1.0,1.0); //white
			//Ogre::ColourValue C(0.0,0.0,0.0); //black

			manual->begin("BaseWhiteNoLighting", Ogre::RenderOperation::OT_LINE_LIST);
			int ne = (int) evert_idx_1.size();
			for (int edge_no=0; edge_no < ne; edge_no++) { //loop over edges
				float pos[3];
				
				vrml.getVertex(part_no, evert_idx_1[edge_no], pos[0], pos[1], pos[2]);
				if (flip_dim > 0)
					pos[flip_dim] = -pos[flip_dim];
				manual->position(pos[0],pos[1],pos[2]);
				manual->colour(C);

				vrml.getVertex(part_no, evert_idx_2[edge_no], pos[0], pos[1], pos[2]);
				if (flip_dim > 0)
					pos[flip_dim] = -pos[flip_dim];
				manual->position(pos[0],pos[1],pos[2]);
				manual->colour(C);
			}
			manual->end();
		}
	}

	node_to_attach_to->attachObject(manual);
}

//
void WmrAnimation::addEntityBox(Ogre::SceneNode* node_to_attach_to, const Ogre::Vector3 L, 
	const Ogre::Quaternion* rot_off, const Ogre::Vector3* trans_off, const Ogre::ColourValue* C) {

	const int nv=8; //number of vertices
	const int nf=6; //number of faces

	Ogre::Vector3 vertices[nv];
	vertices[0] = Ogre::Vector3(-0.5*L[0], -0.5*L[1],  0.5*L[2]);
	vertices[1] = Ogre::Vector3( 0.5*L[0], -0.5*L[1],  0.5*L[2]);
	vertices[2] = Ogre::Vector3( 0.5*L[0],  0.5*L[1],  0.5*L[2]);
	vertices[3] = Ogre::Vector3(-0.5*L[0],  0.5*L[1],  0.5*L[2]);
	vertices[4] = Ogre::Vector3(-0.5*L[0], -0.5*L[1], -0.5*L[2]);
	vertices[5] = Ogre::Vector3( 0.5*L[0], -0.5*L[1], -0.5*L[2]);
	vertices[6] = Ogre::Vector3( 0.5*L[0],  0.5*L[1], -0.5*L[2]);
	vertices[7] = Ogre::Vector3(-0.5*L[0],  0.5*L[1], -0.5*L[2]);

	int indices[nf][4] = {{0,1,2,3}, 
						{5,4,7,6},
						{3,2,6,7},
						{2,1,5,6},
						{1,0,4,5},
						{0,3,7,4}};

	Ogre::Vector3 normals[nf];
	normals[0] = Ogre::Vector3(0,0,1);
	normals[1] = Ogre::Vector3(0,0,-1);
	normals[2] = Ogre::Vector3(0,1,0);
	normals[3] = Ogre::Vector3(1,0,0);
	normals[4] = Ogre::Vector3(0,-1,0);
	normals[5] = Ogre::Vector3(-1,0,0);

	if (rot_off != 0) {
		for (int i=0; i<nv; i++)
			vertices[i] = *rot_off * vertices[i];
		for (int i=0; i<nf; i++)
			normals[i] = *rot_off * normals[i];
	}

	if (trans_off != 0) {
		for (int i=0; i<nv; i++)
			vertices[i] = *trans_off + vertices[i];
	}

	Ogre::ManualObject* manual = mSceneMgr->createManualObject();

	manual->setDynamic(true); //for growing number of vertices

	//vertices with duplicate position, different normal for proper lighting
	manual->estimateVertexCount(nf*6);
	manual->estimateIndexCount(nf*6);

	manual->begin("MyMaterials/FlatVertexColor", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	//add vertices
	for (int face_no=0; face_no < nf; face_no++) { //loop over faces
		for (int vno=0; vno < 4; vno++) { //loop over vertices (in quad)
			manual->position( vertices[indices[face_no][vno]] );
			manual->normal( normals[face_no] );
			if (C != 0)
				manual->colour(*C);
		}
	}

	//add indices, one quad per face
	for (int face_no=0; face_no < nf; face_no++) {
		int vi = face_no*4; //vertex index in the entire list
		manual->quad(vi+0, vi+1, vi+2, vi+3);
	}

	manual->end();

	node_to_attach_to->attachObject(manual);

}

//arrow along z axis (box line + pyramid head)
//length:	total length
//width:	line width
//hlength:	head length
//hwidth:	head width
void WmrAnimation::addEntityArrow(Ogre::SceneNode* node_to_attach_to, 
	const float length, const float width, const float hlength, const float hwidth, 
	const Ogre::Quaternion* rot_off, const Ogre::Vector3* trans_off, const Ogre::ColourValue* C) {


	//make box for line
	Ogre::Vector3 L(width,width,length-hlength);

	Ogre::Vector3 trans_box(0,0,0.5*(length-hlength));
	if (rot_off != 0)
		trans_box = *rot_off * trans_box;
	if (trans_off != 0)
		trans_box = trans_box + *trans_off;

	addEntityBox(node_to_attach_to, L, rot_off, &trans_box, C);

	//make pyramid for head
	//TODO, move this to subfunction?
	
	const int nv = 5;
	const int nf = 6; //6 triangular faces (2 in bottom quad)
	
	Ogre::Vector3 vertices[nv];
	vertices[0] = Ogre::Vector3(-0.5*hwidth, -0.5*hwidth,  length-hlength);
	vertices[1] = Ogre::Vector3( 0.5*hwidth, -0.5*hwidth,  length-hlength);
	vertices[2] = Ogre::Vector3( 0.5*hwidth,  0.5*hwidth,  length-hlength);
	vertices[3] = Ogre::Vector3(-0.5*hwidth,  0.5*hwidth,  length-hlength);
	vertices[4] = Ogre::Vector3(0.0, 0.0, length);

	//anti-clockwise order!
	//pad with -1 for triangular faces
	int indices[nf][3] = {{0,3,2}, //quad
						{0,2,1}, //quad
						{0,1,4},
						{1,2,4},
						{2,3,4},
						{3,0,4}};

	float ang = atan((hwidth/2)/hlength);
	float sa = sin(ang);
	float ca = cos(ang);

	Ogre::Vector3 normals[nf];
	normals[0] = Ogre::Vector3(0,0,-1.0); //quad
	normals[1] = Ogre::Vector3(0,0,-1.0); //quad
	normals[2] = Ogre::Vector3(0,-ca,sa);
	normals[3] = Ogre::Vector3(ca,0,sa);
	normals[4] = Ogre::Vector3(0,ca,sa);
	normals[5] = Ogre::Vector3(-ca,0,sa);

	//rotation, translation
	if (rot_off != 0) {
		for (int i=0; i<nv; i++)
			vertices[i] = *rot_off * vertices[i];
		for (int i=0; i<nf; i++)
			normals[i] = *rot_off * normals[i];
	}

	if (trans_off != 0) {
		for (int i=0; i<nv; i++)
			vertices[i] = *trans_off + vertices[i];
	}

	Ogre::ManualObject* manual = mSceneMgr->createManualObject();

	manual->setDynamic(true); //for growing number of vertices

	manual->estimateVertexCount(nf*3);
	manual->estimateIndexCount(nf*3);

	manual->begin("MyMaterials/FlatVertexColor", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	
	//add vertices
	for (int face_no=0; face_no<nf; face_no++) { //loop over faces
		for (int vno=0; vno<3; vno++) { //loop over vertices (in triangle)
			manual->position(vertices[indices[face_no][vno]]);
			manual->normal(normals[face_no]);
			if (C != 0)
				manual->colour(*C);
		}
	}

	//add indices, one quad per face
	for (int face_no=0; face_no<nf; face_no++) {
		int vi = face_no*3; //vertex index in the entire list
		manual->triangle(vi+0, vi+1, vi+2);
	}
	manual->end();

	node_to_attach_to->attachObject(manual);
	
}

void WmrAnimation::addEntitySurface(Surface* surf) {
	TriMeshSurf* trimesh_surf = dynamic_cast<TriMeshSurf*>(surf);
	if (trimesh_surf != 0) {
		addEntityTriMeshSurf(trimesh_surf);
		return;
	}

	GridSurf* grid_surf = dynamic_cast<GridSurf*>(surf);
	if (grid_surf != 0) {
		addEntityGridSurf(grid_surf);
		return;
	}
}

void WmrAnimation::addEntityTriMeshSurf(TriMeshSurf* surf) {
	
	//get from TriMeshSurf object
	const int nt = surf->get_nt(); //number of triangles
	const Vec3* vertices = surf->get_vertices();
	const int* indices = surf->get_indices();
	const Real* pec = surf->get_pec();

	Ogre::ManualObject* manual = mSceneMgr->createManualObject();

	manual->setDynamic(true); //for growing number of vertices
	//before begin or after?
	manual->estimateVertexCount(3*nt);
	manual->estimateIndexCount(3*nt);

	manual->begin("MyMaterials/FlatVertexColor", Ogre::RenderOperation::OT_TRIANGLE_LIST);
	Ogre::ColourValue C(0.0,1.0,1.0);

	for (int tno=0; tno<nt; tno++) { //loop over triangles
		for (int vno=0; vno<3; vno++) { //loop over vertices (in triangle)
			const Real* vert = vertices[indices[S2I(vno,tno,3)]]; //pointer to vertex
			const Real* pec_ = pec + S2I(0,tno,4); //pointer to normal

			manual->position(vert[0], vert[1], vert[2]);
			manual->normal(pec_[0],pec_[1],pec_[2]);
			manual->colour(C);
		}
	}

	for (int tno=0; tno<nt; tno++) {
		int vi = tno*3; //vertex index in the entire list
		manual->triangle(vi+0, vi+1, vi+2);
		manual->triangle(vi+0, vi+2, vi+1); //clockwise order to view from bottom
	}

	manual->end();
	mSceneMgr->getRootSceneNode()->attachObject(manual);

}

void WmrAnimation::addEntityGridSurf(GridSurf* surf) {

	//get from GridSurf object
	const float lowx = surf->get_lowerlimx();
	const float lowy = surf->get_lowerlimy();
	const int nx = surf->get_nx();
	const int ny = surf->get_ny();
	const float dx = (float) surf->get_dx();
	const float dy = (float) surf->get_dy();

	Ogre::ManualObject* manual = mSceneMgr->createManualObject();

	manual->setDynamic(true); //for growing number of vertices
	//before begin or after?
	manual->estimateVertexCount(nx*ny);
	manual->estimateIndexCount(6*(nx-1)*(ny-1));
	//materials in media/materials/scripts directory of sdk

	//manual->begin("MyMaterials/FlatVertexColor", Ogre::RenderOperation::OT_TRIANGLE_LIST);
	//Ogre::ColourValue C(0.0,1.0,1.0);
	
	manual->begin("MyMaterials/Terrain", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	Vec3 pos; //position
	Vec3 nvec; //normal vector
	const Real* Z = surf->get_Z();

	for (int j=0; j<ny; j++) {
		for (int i=0; i<nx; i++) {

			pos[0] = i*dx+lowx;
			pos[1] = j*dy+lowy;

			pos[2] = Z[i + j*nx];
			surf->surfaceNormal(pos, true, nvec);

			manual->position(pos[0], pos[1], pos[2]);

			manual->normal(nvec[0], nvec[1], nvec[2]);
			manual->textureCoord(float(i)/float(nx-1), float(j)/float(ny-1));
			//manual->colour(C);
		}
	}

	for (int i=0; i<nx-1; i++) {
		for (int j=0; j<ny-1; j++) {
			manual->quad(i+j*nx, (i+1)+j*nx, (i+1)+(j+1)*nx, i+(j+1)*nx); //counter-clockwise
		}
	}

	manual->end();

	//change scale here or in .material script
	manual->getSection(0)->getMaterial()->getTechnique(0)->getPass(0)->getTextureUnitState(0)->setTextureScale(.1,.1);

	mSceneMgr->getRootSceneNode()->attachObject(manual);
		
}

void WmrAnimation::addLine( const int num_chains ) {
	//use billboard chain to draw thick lines
	assert(nl<MAXNL);

	lines[nl] = mSceneMgr->createBillboardChain();
	
	lines[nl]->setNumberOfChains(num_chains);
	lines[nl]->setMaxChainElements(1e3);
	lines[nl]->setMaterialName("MyMaterials/BaseWhiteNoLighting"); //cull_hardware none to be visible from both sides

	//Lines[nl]->setFaceCamera(true); //not working! have to do manually in render loop
	mSceneMgr->getRootSceneNode()->attachObject(lines[nl]);

	//size_t nce = Lines[nl]->getNumChainElements(0); //DEBUGGING, why 1 not 0?

	nl++;

	
}

void WmrAnimation::addChainElement( const int lineno, const int chainno, const Ogre::Vector3 pt, const float thickness, Ogre::ColourValue* C) {
	
	Ogre::BillboardChain::Element elem;
	elem.position = pt;

	size_t nce = lines[lineno]->getNumChainElements(chainno);
	if (nce>1) { //why 1 not 0?
		//make appearance match first element
		elem.width = lines[lineno]->getChainElement(chainno, 0).width;
		elem.colour = lines[lineno]->getChainElement(chainno, 0).colour;
	} else { //first element
		elem.width = thickness;
		if (C != 0)
			elem.colour = *C;
	}

	lines[lineno]->addChainElement( chainno, elem);
}

void WmrAnimation::updateNodesLines( const int nf, const HomogeneousTransform HT_parent[], const int num_contacts, const ContactGeom* contacts ) {
	//note assumptions made!

	for (int fi=0; fi < nf; fi++) //loop over nodes
		updateNodeTransform(fi,HT_parent[fi]);

	//lines
	//body frame
	Ogre::ColourValue C(0.0,0.0,1.0);
	float thk = .01f;
	Ogre::Vector3 pt = nodes[0]->_getDerivedPosition(); //assumes Nodes[0] for body frame
	addChainElement(0,0,pt,thk,&C); //assumes Lines[0] for body frame

	
	const WheelContactGeom* wcontacts = dynamic_cast<const WheelContactGeom*>(contacts);

	if (wcontacts != 0) {
		const int nw = num_contacts;

		//contact points
		C = Ogre::ColourValue(1.0,0.0,0.0);
		thk = .01f;

		for (int wno=0; wno < nw; wno++) {
			int node_ind = nf + wno; //node index for contact frame
			int line_ind = 1;

			//update transform
			updateNodeTransform(node_ind,wcontacts[wno].HT_wheel);
			nodes[node_ind]->setVisible(wcontacts[wno].incontact); //only visible if in contact

			//update lines
			Real dz = wcontacts[wno].dz;

			pt = nodes[nf+wno]->_getDerivedPosition();
			pt = pt - (dz-thk/2) * nodes[node_ind]->_getDerivedOrientation().zAxis(); //move to surface
			addChainElement(line_ind,wno,pt,thk,&C);
		}

		return;
	} 

	//TODO
	const TrackContactGeom* tcontacts = dynamic_cast<const TrackContactGeom*>(contacts);

	if (tcontacts != 0) {
		const int nt = num_contacts;

		//set visibility of contact frame arrows
		int node_idx = nf-1;
		for (int tno=0; tno < nt; tno++) {
			node_idx++; //for track node
			int np = tcontacts[tno].get_np();
			for (int pno=0; pno < np; pno++) {
				node_idx++;
				nodes[node_idx]->setVisible(tcontacts[tno].incontact[pno]);
			}
		}
	}
}

void WmrAnimation::go(void) {

	//copied from http://www.ogre3d.org/tikiwiki/tiki-index.php?page=Basic+Tutorial+6#Render_Loop_-_take_one
	// Infinite loop, until broken out of by frame listeners
    while( 1 )
    {
        //Pump messages in all registered RenderWindow windows
        Ogre::WindowEventUtilities::messagePump();


        if (!mRoot->renderOneFrame())
            break;
		
    }
}

void WmrAnimation::start(void) {
#ifdef _DEBUG
	mResourcesCfg = "resources_d.cfg";
	mPluginsCfg = "plugins_d.cfg";
#else
	mResourcesCfg = "resources.cfg";
	mPluginsCfg = "plugins.cfg";
#endif
	
	if (!setup())
        return;
}

bool WmrAnimation::updateRender(void) {

	for (int i=0; i<nl; i++) 
		lines[i]->setFaceCamera(false, mCamera->getDirection()); 

	//Pump messages in all registered RenderWindow windows
	Ogre::WindowEventUtilities::messagePump();

	return mRoot->renderOneFrame();

	
}
