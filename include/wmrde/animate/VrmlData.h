//VrmlData.h
//This class parses vertex, index, and color information from VRML files (.wrl) for use in animation

#ifndef _VRML_DATA_H_
#define _VRML_DATA_H_

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>


template<typename Type> //Type is float or int
void stringToNumVector( std::string& tmp, std::vector<Type>& out ) { 
	//remove commas, line breaks such that numbers are delimited only by spaces
	tmp.erase(std::remove(tmp.begin(), tmp.end(), ','), tmp.end());
	tmp.erase(std::remove(tmp.begin(), tmp.end(), '\n'), tmp.end());
	
	std::stringstream ss(tmp);
	while (!ss.fail()) {
		Type f;
		ss >> f;
		if (!ss.fail()) 
			out.push_back(f);
	}
}

template <typename T> //float or double, TODO a better way?
class VrmlData {
public:
	static const int MAX_NUM_PARTS = 20;
private:
	std::vector<T> coord[MAX_NUM_PARTS];
	std::vector<int> coord_idx[MAX_NUM_PARTS];
	std::vector<T> normal[MAX_NUM_PARTS];
	std::vector<int> normal_idx[MAX_NUM_PARTS];
	std::vector<T> color[MAX_NUM_PARTS];
	std::vector<int> color_idx[MAX_NUM_PARTS];
	int nparts_;
public:
	VrmlData() { }

	int readfile(const std::string filename) {
		// declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
	
		std::ifstream file ( filename.c_str() );

		int part_no=-1;

		std::string tmp;
		if ( file.is_open() ) {
			while ( file.good() ) {
				std::getline(file, tmp, '\n');

				if ( tmp.find("color [") < tmp.length() ) {
					if (part_no+1 >= MAX_NUM_PARTS)
						break;
					part_no++;
					std::getline(file, tmp, ']');
					stringToNumVector( tmp, color[part_no] );
				}

				if ( tmp.find("point [") < tmp.length() ) {
					std::getline(file, tmp, ']');
					stringToNumVector( tmp, coord[part_no] );
				}

				if ( tmp.find("vector [") < tmp.length() ) {
					std::getline(file, tmp, ']');
					stringToNumVector( tmp, normal[part_no] );
				}

				if ( tmp.find("coordIndex [") < tmp.length() ) {
					std::getline(file, tmp, ']');
					stringToNumVector( tmp, coord_idx[part_no] );
				}

				if ( tmp.find("colorIndex [") < tmp.length() ) {
					std::getline(file, tmp, ']');
					stringToNumVector( tmp, color_idx[part_no] );
				}

				if ( tmp.find("normalIndex [") < tmp.length() ) {
					std::getline(file, tmp, ']');
					stringToNumVector( tmp, normal_idx[part_no] );
				}
			}
			file.close();
			nparts_ = part_no+1;
			return 1;
		} else {
			std::cout << "could not open .wrl file" << std::endl;
			return 0;
		}
	}

	int nparts() { return nparts_; }
	int ncoords(int part_no) { return (int) coord[part_no].size() / 3; }
	int ntriangles(int part_no) { return (int) coord_idx[part_no].size() / 4; } 

	int getTriangleVertexIndex(int part_no, int tri_no, int vert_no) { return coord_idx[part_no][tri_no*4 + vert_no]; }
	int getTriangleColorIndex(int part_no, int tri_no) {
		int idx=0; //in case no color_idx
		if (color_idx[part_no].size() > 0)
			idx = color_idx[part_no][tri_no];
		return idx;
	}
	int getTriangleNormalIndex(int part_no, int tri_no) { return normal_idx[part_no][tri_no*4+0]; }

	void getVertex(int part_no, int idx, T& x, T& y, T & z) {
		x = coord[part_no][idx*3 + 0];
		y = coord[part_no][idx*3 + 1];
		z = coord[part_no][idx*3 + 2];
	}
	void getColor(int part_no, int idx, T& R, T& G, T& B) {
		R=color[part_no][idx*3 + 0];
		G=color[part_no][idx*3 + 1];
		B=color[part_no][idx*3 + 2];
	}
	void getNormal(int part_no, int idx, T& nx, T& ny, T& nz) {
		nx = normal[part_no][idx*3 + 0];
		ny = normal[part_no][idx*3 + 1];
		nz = normal[part_no][idx*3 + 2];
	}
	
	void getTriangleVertex(int part_no, int tri_no, int vert_no, T& x, T& y, T& z) {
		getVertex( part_no, getTriangleVertexIndex(part_no,tri_no,vert_no), x,y,z);
	}
	void getTriangleColor(int part_no, int tri_no, T& R, T& G, T& B) {
		getColor( part_no, getTriangleColorIndex(part_no, tri_no), R,G,B);
	}
	void getTriangleNormal(int part_no, int tri_no, T& nx, T& ny, T& nz) {
		getNormal( part_no, getTriangleNormalIndex(part_no,tri_no), nx,ny,nz);
	}

};

#endif


