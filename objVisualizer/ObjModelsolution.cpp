/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "ObjModel.hpp"



#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;


ObjModel::ObjModel()
{
	_numVertices = 0;
	_numTriangles = 0;
}

/**
 * Calculate the normal of a triangular face defined by three points
 * 
 * @param[in] coord1 the first vertex
 * @param[in] coord2 the second vertex
 * @param[in] coord3 the third vertex
 * @param[out] norm the normal
 */
void ObjModel::computeNormal( const float coord1[3], const float coord2[3], const float coord3[3], float norm[3]  ) const
{
	//*********************************************
	// Compute the normal vector of the 3 input vertices:
	// Considering the first vertex as "reference", compute the 
	// vectors connecting the reference and the other 2 vertices.
	// The normal is given by the cross product of these two vertices
	//
	// Remember to normalize the final vector
	//*********************************************
	float va[3], vb[3], val;
	
	//**************************************************
	// Compute the vector connecting the "reference" with the second vertex
	//**************************************************
	va[0] = coord1[0] - coord2[0];
	va[1] = coord1[1] - coord2[1];
	va[2] = coord1[2] - coord2[2];

	//**************************************************
	// Compute the vector connecting the "reference" with the third vertex
	//**************************************************
	vb[0] = coord1[0] - coord3[0];
	vb[1] = coord1[1] - coord3[1];
	vb[2] = coord1[2] - coord3[2];

	//**************************************************
	// Compute the normal as the cross product of the latters
	//**************************************************
	norm[0] = va[1] * vb[2] - vb[1] * va[2];
	norm[1] = vb[0] * va[2] - va[0] * vb[2];
	norm[2] = va[0] * vb[1] - vb[0] * va[1];

	//**************************************************
	// remember to normalize the vector
	//**************************************************
	val = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );

	norm[0] = norm[0]/val;
	norm[1] = norm[1]/val;
	norm[2] = norm[2]/val;
 
}

void ObjModel::computeNormal( const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm  ) const
{
	norm = (v1-v2).cross( v1-v3 );

	norm.normalize();
}

/**
 * Perform a first scan of the file in order to get the number of vertices and the number of faces
 * @param[in] filename the OBJ file
 * @param[out] vertexNum the number of vertices found
 * @param[out] triangleNum the number of faces found
 */
void ObjModel::firstScan(char* filename, long &vertexNum, long &triangleNum)
{
	string line;
	ifstream objFile (filename);	
	// If obj file is open, continue
	if (objFile.is_open())													
	{
		vertexNum = 0;
		triangleNum = 0;
		
		// Start reading file data
		while (! objFile.eof() )											
		{	
			//**************************************************
			// get a line of the file (use getline )
			//**************************************************
			getline (objFile,line);											
 
			//**************************************************
			// If the first character is a 'v'...
			//**************************************************
			if (line.c_str()[0] == 'v')										
			{
				//**************************************************
				// Increment the number of vertices
				//**************************************************
				++vertexNum;
			}
			//**************************************************
			// If the first character is a 'f' ...
			//**************************************************
			if (line.c_str()[0] == 'f')										
			{
				//**************************************************
				// Increment the number of triangles
				//**************************************************
				++triangleNum;
			}
		}
		// Close OBJ file
		objFile.close();														
	}
	else 
	{
		cout << "Unable to open file";								
	}
}



/**
 * @brief 
 * @param filename
 * @return 
 */
int ObjModel::load(char* filename)
{
	// Count the number of vertices and the number of triangles
	firstScan(filename, _numVertices, _numTriangles);
	
cout << "Found object with "<< _numVertices << " vertices and "  << _numTriangles << " faces" << endl;

	string line;
	ifstream objFile (filename);	
	
	// If obj file is open, continue
	if (objFile.is_open())													
	{
		
		vector<int> vOcc;
 
		//**************************************************
		// Allocate memory for the vertices: use malloc to allocate the memory
		// How many floats do you need overall...?
		//**************************************************
		_vertices = (float*) malloc ( _numVertices*COORD_PER_VERTEX*sizeof(float) );	
		
		//**************************************************
		// Allocate memory for the triangles: use malloc to allocate the memory
		// How many floats do you need overall...? count how many floats you need
		// for each vertices, hence how many floats for each triangle...
		//**************************************************
		_triangles = (float*) malloc(_numTriangles*TOTAL_FLOATS_IN_TRIANGLE*2*sizeof(float));			
		//**************************************************
		// Allocate memory for the triangles: use normals to allocate the memory
		// How many floats do you need overall...?
		//**************************************************		
		_normals  = (float*) malloc(_numTriangles*TOTAL_FLOATS_IN_TRIANGLE*sizeof(float));					
		
		
		// This index is used to run through the triangle array
		int triangle_index = 0;	
		// This index is used to run through the normal array
		int normal_index = 0;		
		// This index is used to run through the vertex array
		int vertex_index = 0;
 
		// Start reading file data
		while (! objFile.eof() )											
		{		
			//**************************************************
			// Get a line from file
			//**************************************************
			getline (objFile,line);											
 
			//**************************************************
			// If the first character is a 'v'...
			//**************************************************
			if (line.c_str()[0] == 'v')										
			{
//			    cout << line.c_str()[0] << endl;
				
				// Set first character to ' '. This will allow us to use sscanf
				line[0] = ' ';												
 
				//**************************************************
				// Read 3 floats from the line:  X Y Z and store them in the corresponding place in _vertices
				// In order to read the floats in one shot from string you can use sscanf
				//**************************************************
				sscanf(line.c_str(),"%f %f %f ",							
					&_vertices[vertex_index],
					&_vertices[vertex_index+1], 
					&_vertices[vertex_index+2]);

				_v.push_back(point3d(_vertices[vertex_index],_vertices[vertex_index+1],_vertices[vertex_index+2]));
				vOcc.push_back(0);
				_nv.push_back(vec3d());
				
				//**************************************************
				// This is for the 2nd part of the exercise: update the bounding box
				// For the very first vertex read, initialize the bb accordingly
				//**************************************************
				if( vertex_index == 0)
				{
					//**************************************************
					// Case of the very first vertex read
					//**************************************************
					_bb.Xmax = _vertices[vertex_index];
					_bb.Xmin = _vertices[vertex_index];
					_bb.Ymax = _vertices[vertex_index+1];
					_bb.Ymin = _vertices[vertex_index+1];
					_bb.Zmax = _vertices[vertex_index+2];
					_bb.Zmin = _vertices[vertex_index+2];
				}
				else
				{
					//**************************************************
					// otherwise check the vertex against the bounding box and in case update it
					//**************************************************
					_bb.Xmax = ( _bb.Xmax < _vertices[vertex_index] )	? (_vertices[vertex_index]) : (_bb.Xmax) ;
					_bb.Xmin = ( _bb.Xmin > _vertices[vertex_index] )	? (_vertices[vertex_index]) : (_bb.Xmin) ;		
					_bb.Ymax = ( _bb.Ymax < _vertices[vertex_index+1] )	? (_vertices[vertex_index+1]) : (_bb.Ymax) ;
					_bb.Ymin = ( _bb.Ymin > _vertices[vertex_index+1] )	? (_vertices[vertex_index+1]) : (_bb.Ymin) ;		
					_bb.Zmax = ( _bb.Zmax < _vertices[vertex_index+2] )	? (_vertices[vertex_index+2]) : (_bb.Zmax) ;
					_bb.Zmin = ( _bb.Zmin > _vertices[vertex_index+2] )	? (_vertices[vertex_index+2]) : (_bb.Zmin) ;		
				}
					
				// update the vertex
				vertex_index += COORD_PER_VERTEX;		
				
				// just a security check...
				assert((vertex_index <= _numVertices*3));
			}
			//**************************************************
			// If the first character is a 'f'...
			//**************************************************
			if (line.c_str()[0] == 'f')										
			{
//				cout << line.c_str()[0] << endl;
				
				// Set first character to ' '. This will allow us to use sscanf
				line[0] = ' ';												
 
				// this contains temporary the indices of the vertices
				int vertexIdx[3] = { 0, 0, 0 };
				
				//**************************************************
				// Read 3 integers from the line:  idx1 idx2 idx3 and store them in the corresponding place in vertexIdx
				// In order to read the 3 integers in one shot from string you can use sscanf
				//**************************************************
				sscanf(line.c_str(),"%i %i %i",								// Read integers from the line:  f 1 2 3
					&vertexIdx[0],										// First point of our triangle. This is an 
					&vertexIdx[1],										// pointer to our vertexBuffer list
					&vertexIdx[2] );										// each point represents an X,Y,Z.
 

				
				//**************************************************
				// correct the indices: OBJ starts counting from 1, in C the arrays starts at 0...
				//**************************************************
				vertexIdx[0] -= 1;
				vertexIdx[1] -= 1;		
				vertexIdx[2] -= 1;										

				 //cout << vertexNumber[0] << " " << vertexNumber[1] << " " << vertexNumber[2] << " "  << endl;
				_indices.push_back(triangleIndex(vertexIdx[0],vertexIdx[1],vertexIdx[2]));
				vOcc[vertexIdx[0]]++;
				vOcc[vertexIdx[1]]++;
				vOcc[vertexIdx[2]]++;
				
				//just a security check
				assert((vertexIdx[0] >= 0 )&&(vertexIdx[0] <= _numVertices));
				assert((vertexIdx[1] >= 0 )&&(vertexIdx[1] <= _numVertices));
				assert((vertexIdx[2] >= 0 )&&(vertexIdx[2] <= _numVertices));
								
				//*********************************************************************
				//  fill the _triangles array with the 3 vertices
				// tCounter gives you the starting position of each vertex in _triangles
				//*********************************************************************
 				int tCounter = 0;
				for (int i = 0; i < VERTICES_PER_TRIANGLE; i++)					
				{
					_triangles[triangle_index + tCounter   ] = _vertices[COORD_PER_VERTEX*vertexIdx[i] ];
					_triangles[triangle_index + tCounter +1 ] = _vertices[COORD_PER_VERTEX*vertexIdx[i]+1 ];
					_triangles[triangle_index + tCounter +2 ] = _vertices[COORD_PER_VERTEX*vertexIdx[i]+2 ];
					
					tCounter += VERTICES_PER_TRIANGLE;
				}
 
				//*********************************************************************
				//  Calculate the normal of the triangles, it will be the same for each vertex
				//*********************************************************************
				float norm[3];
				
				//*********************************************************************
				//  compute the normal for the 3 vertices we just added
				//*********************************************************************
				computeNormal( &_vertices[COORD_PER_VERTEX*vertexIdx[0] ], 
									&_vertices[COORD_PER_VERTEX*vertexIdx[1] ], 
									&_vertices[COORD_PER_VERTEX*vertexIdx[2] ], 
									norm );

				cout << "n " << norm[0] << " " << norm[1] << " " << norm[2] << endl;
 
				
				//*********************************************************************
				//  fill the _normals array with 3 copy of the normal
				// tCounter gives you the starting position of each normal in _normals
				//*********************************************************************
				tCounter = 0;
				for (int i = 0; i < VERTICES_PER_TRIANGLE; i++)
				{
					_normals[normal_index + tCounter ] = norm[0];
					_normals[normal_index + tCounter +1] = norm[1];
					_normals[normal_index + tCounter +2] = norm[2];
					tCounter += VERTICES_PER_TRIANGLE;
					_nt.push_back(vec3d(norm[0],norm[1],norm[2]));
				}
 

				PRINTVAR(angleAtVertex(_v[vertexIdx[0]], _v[vertexIdx[1]], _v[vertexIdx[2]]));
				PRINTVAR(angleAtVertex(_v[vertexIdx[1]], _v[vertexIdx[0]], _v[vertexIdx[2]]));
				PRINTVAR(angleAtVertex(_v[vertexIdx[2]], _v[vertexIdx[0]], _v[vertexIdx[1]]));
				_nv[vertexIdx[0]] += (vec3d(norm) * angleAtVertex(_v[vertexIdx[0]], _v[vertexIdx[1]], _v[vertexIdx[2]]));
				_nv[vertexIdx[1]] += (vec3d(norm) * angleAtVertex(_v[vertexIdx[1]], _v[vertexIdx[0]], _v[vertexIdx[2]]));
				_nv[vertexIdx[2]] += (vec3d(norm) * angleAtVertex(_v[vertexIdx[2]], _v[vertexIdx[0]], _v[vertexIdx[1]]));
				
				// update the indices
				triangle_index += TOTAL_FLOATS_IN_TRIANGLE;
				normal_index += TOTAL_FLOATS_IN_TRIANGLE;
	
				
				
				// just a security check
				assert((triangle_index <= _numTriangles*TOTAL_FLOATS_IN_TRIANGLE));		
				assert((normal_index <= _numTriangles*TOTAL_FLOATS_IN_TRIANGLE));

			}	
		}
		cout << "TotalConnectedTriangles "<< triangle_index << endl;

		cout << "Arrays:\n\tindices "<< _indices.size() << "\t_v " << _v.size() << "\t_nt " << _nt.size() << "\t_nv " << _nv.size() << endl;
		PRINTVAR( _indices );
		PRINTVAR( _v );
		PRINTVAR( _nt );
		PRINTVAR( _nv );
		for(int i=0; i<_nv.size(); _nv[i++].normalize());
		PRINTVAR( _nv );
		
		assert( _numTriangles == triangle_index/TOTAL_FLOATS_IN_TRIANGLE );
		assert( _numVertices == vertex_index/COORD_PER_VERTEX );
		
		
		// Close OBJ file
		objFile.close();

	}
	else 
	{
		cout << "Unable to open file";								
	}
	

	cout << "Object loaded with "<< _numVertices << " vertices and "  << _numTriangles << " faces" << endl;
	cout << "Bounding box : Xmax=" << _bb.Xmax << "  Xmin=" << _bb.Xmin << "  Ymax=" << _bb.Ymax << "  Ymin=" << _bb.Ymin << "  Zmax=" << _bb.Zmax << "  Zmin=" << _bb.Zmin << endl;
	return 0;
}


/**
 * Computes the angle at vertex baseV formed by the edges connecting it with the
 * vertices v1 and v2 respectively, ie the baseV-v1 and baseV-v2 edges
 * @brief
 * @param baseV the vertex at which to compute the angle
 * @param v1 the other vertex of the first edge baseV-v1
 * @param v2 the other vertex of the second edge baseV-v2
 * @return the angle in radiants
 */
float ObjModel::angleAtVertex( const point3d& baseV, const point3d& v1, const point3d& v2 ) const
{
	vec3d e1 = baseV-v1;
	vec3d e2 = baseV-v2;
//	PRINTVAR(e1);
//	PRINTVAR(e2);
//	PRINTVAR(e1.norm());
//	PRINTVAR(e2.norm());
//	PRINTVAR((e1.norm()*e1.norm()));
//	PRINTVAR((e1).dot(e2) / (e1.norm()*e1.norm()));
	return ( acos( (e1).dot(e2) / (e1.norm()*e2.norm()) ) );
}


void ObjModel::subdivision()
{
    // for each triangle
        // for each edge

}

GLushort ObjModel::getNewVertex( const edge &e, vector<point3d> &vertList, EdgeList &newVertList ) const
{
    // if the egde is not contained in the new vertex list
    if( newVertList.contains( e ) )
    {
        // generate new index (vertex.size)
        GLushort idxnew = vertList.size();
        // add the edge and index to the map
        newVertList.add(e , idxnew );
        // generate new vertex
        point3d nvert = (vertList[e.first] + vertList[e.second])*0.5;
        // append it to the list of vertex
        vertList.push_back( nvert );
        // normals?

    }
    // else
    {
        // get the index of the vertex
    }
}

void ObjModel::release()
{
	free(this->_triangles);
	free(this->_normals);
	free(this->_vertices);
}
 
void ObjModel::draw() const
{
	glShadeModel( GL_SMOOTH );
	
	//****************************************
	// Enable vertex arrays
	//****************************************
	glEnableClientState(GL_VERTEX_ARRAY);
	
	//****************************************
	// Enable normal arrays
	//****************************************
	glEnableClientState(GL_NORMAL_ARRAY);
	
	//****************************************
	// Vertex Pointer to triangle array
	//****************************************
	glVertexPointer(3,GL_FLOAT,	0,_triangles);
	
	//****************************************
	// Normal pointer to normal array
	//****************************************
	glNormalPointer(GL_FLOAT, 0, _normals);	
	
	//****************************************
	// Draw the triangles
	//****************************************
	glDrawArrays(GL_TRIANGLES, 0, _numTriangles*VERTICES_PER_TRIANGLE);	

	//****************************************
	// Disable vertex arrays	
	//****************************************
	glDisableClientState(GL_VERTEX_ARRAY);
	
	//****************************************
	// Disable normal arrays
	//****************************************
	glDisableClientState(GL_NORMAL_ARRAY);					
	
	
//	for(int i=0; i<_numTriangles; ++i)
//	{
//		glBegin(GL_TRIANGLES);
//			glNormal3fv(&_normals[i*TOTAL_FLOATS_IN_TRIANGLE]);
//			glVertex3fv(&_triangles[i*TOTAL_FLOATS_IN_TRIANGLE]);
//
//			glVertex3fv(&_triangles[i*TOTAL_FLOATS_IN_TRIANGLE+1*COORD_PER_VERTEX]);
//
//			glVertex3fv(&_triangles[i*TOTAL_FLOATS_IN_TRIANGLE+2*COORD_PER_VERTEX]);
//		glEnd();
//	}

}

void ObjModel::flatDraw() const
{
	glShadeModel( GL_SMOOTH );

	// for each triangle draw the vertices and the normals
	for(int i=0; i<_indices.size(); i++)
	{
		glBegin(GL_TRIANGLES);
			//compute the normal of the triangle
			vec3d n;
			computeNormal( _v[_indices[i].v1], _v[(int)_indices[i].v2], _v[(int)_indices[i].v3], n);
			glNormal3fv((float*)&n);

			glVertex3fv((float*)&_v[_indices[i].v1]);

			glVertex3fv((float*)&_v[_indices[i].v2]);

			glVertex3fv((float*)&_v[_indices[i].v3]);

		glEnd();
	}

}

void ObjModel::wireframetDraw() const
{

	glDisable(GL_LIGHTING);
	glColor3f(0,0,0);
	glLineWidth(2);
	for(int i=0; i<_indices.size(); i++)
	{
		glBegin(GL_LINE_LOOP);
			glVertex3fv((float*)&_v[_indices[i].v1]);

			glVertex3fv((float*)&_v[_indices[i].v2]);

			glVertex3fv((float*)&_v[_indices[i].v3]);
		glEnd();
	}
	glEnable(GL_LIGHTING);

}

void ObjModel::indexDraw() const
{
	glShadeModel( GL_SMOOTH );

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glNormalPointer(GL_FLOAT, 0, (float*)&_nv[0]);
	glVertexPointer(COORD_PER_VERTEX, GL_FLOAT, 0, (float*)&_v[0]);

	glDrawElements(GL_TRIANGLES, _numTriangles*VERTICES_PER_TRIANGLE, GL_UNSIGNED_SHORT, (GLushort*)&_indices[0]);


	glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);
}


/**
 * It scales the model to unitary size by translating it to the origin and
 * scaling it to fit in a unit cube around the origin.
 * 
 * @return the scale factor used to transform the model
 */
float ObjModel::unitizeModel()
{
	if(_vertices && _triangles)
	{     
		//****************************************
		// calculate model width, height, and 
		// depth using the bounding box
		//****************************************
		float w,h,d;
		w = fabs(_bb.Xmax) + fabs(_bb.Xmin);
		h = fabs(_bb.Ymax) + fabs(_bb.Ymin);
		d = fabs(_bb.Zmax) + fabs(_bb.Zmin);
		
		//****************************************
		// calculate center of the bounding box of the model 
		//****************************************
		float cx = (_bb.Xmax + _bb.Xmin) / 2.0;
		float cy = (_bb.Ymax + _bb.Ymin) / 2.0;
		float cz = (_bb.Zmax + _bb.Zmin) / 2.0;
		
		//****************************************
		// calculate the unitizing scale factor as the 
		// maximum of the 3 dimensions
		//****************************************
		float scale = 2.0 / std::max(std::max(w, h), d);
		
cout << "scale: " << scale << " cx " << cx << " cy " << cy << " cz " << cz << endl;		

		// translate each vertex wrt to the center and then apply the scaling to the coordinate
		for (int i = 0; i < _numVertices; i++) 
		{
			//****************************************
			// translate the vertex
			//****************************************
			_vertices[COORD_PER_VERTEX * i + 0] -= cx;
			_vertices[COORD_PER_VERTEX * i + 1] -= cy;
			_vertices[COORD_PER_VERTEX * i + 2] -= cz;

			_v[i].translate(-cx, -cy, -cz);
			
//		 	cout << vertexBuffer[COORD_PER_VERTEX * i +0] << " " << vertexBuffer[COORD_PER_VERTEX * i +1] << " " << vertexBuffer[COORD_PER_VERTEX * i +2] << " "  << endl;
			
			//****************************************
			// apply the scaling
			//****************************************
			_vertices[COORD_PER_VERTEX * i + 0] *= scale;
			_vertices[COORD_PER_VERTEX * i + 1] *= scale;
			_vertices[COORD_PER_VERTEX * i + 2] *= scale;

			_v[i].scale(scale);
			
//		    cout << vertexBuffer[COORD_PER_VERTEX * i +0] << " " << vertexBuffer[COORD_PER_VERTEX * i +1] << " " << vertexBuffer[COORD_PER_VERTEX * i +2] << " "  << endl;
		}


		// do the same for each vertex of the triangles
		for (int i = 0; i < _numTriangles; ++i) 
		{		
			// for each triangle
			for (int j = 0; j < VERTICES_PER_TRIANGLE; ++j)					
			{
				//****************************************
				// translate the vertex
				//****************************************
				_triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j*COORD_PER_VERTEX   ] -= cx;
				_triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j*COORD_PER_VERTEX +1 ] -= cy;
				_triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j*COORD_PER_VERTEX +2 ] -= cz;
				
//cout << Faces_Triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j +0] << " " << Faces_Triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j +1] << " " << Faces_Triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j +2] << " "  << endl;
				
				//****************************************
				// apply the scaling
				//****************************************
				_triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j*COORD_PER_VERTEX   ] *= scale; 
				_triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j*COORD_PER_VERTEX +1 ] *= scale;
				_triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j*COORD_PER_VERTEX +2 ] *= scale;
//cout << _triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j +0] << " " << _triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j +1] << " " << _triangles[i*TOTAL_FLOATS_IN_TRIANGLE + j +2] << " "  << endl;
			}

		}

		//****************************************
		// update the bounding box, ie translate and scale the 6 coordinates
		//****************************************
		_bb.Xmax = (_bb.Xmax - cx) * scale;
		_bb.Xmin = (_bb.Xmin - cx) * scale;
		_bb.Ymax = (_bb.Ymax - cy) * scale;
		_bb.Ymin = (_bb.Ymin - cy) * scale;
		_bb.Zmax = (_bb.Zmax - cz) * scale;
		_bb.Zmin = (_bb.Zmin - cz) * scale;
		
		cout << "New bounding box : Xmax=" << _bb.Xmax << "  Xmin=" << _bb.Xmin << "  Ymax=" << _bb.Ymax << "  Ymin=" << _bb.Ymin << "  Zmax=" << _bb.Zmax << "  Zmin=" << _bb.Zmin << endl;
		
		return scale;
		
	}
	
	return 0;
}


