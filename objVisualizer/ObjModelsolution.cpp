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
#include <sstream>


using namespace std;


ObjModel::ObjModel()
{
}

/**
 * Calculate the normal of a triangular face defined by three points
 *
 * @param[in] v1 the first vertex
 * @param[in] v2 the second vertex
 * @param[in] cv3 the third vertex
 * @param[out] norm the normal
 */
void ObjModel::computeNormal( const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm  ) const
{
	norm = (v1-v2).cross( v1-v3 );

	norm.normalize();
}

/**
 * @brief 
 * @param filename
 * @return 
 */
int ObjModel::load(char* filename)
{
	
	string line;
	ifstream objFile (filename);	
	
	// If obj file is open, continue
	if (objFile.is_open())													
	{
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
			PRINTVAR(line);
			if ((line.c_str()[0] == 'v') && (line.c_str()[1] == ' '))	// to drop all the vn and vn lines
			{
				PRINTVAR(line);
				// Set first character to ' '. This will allow us to use sscanf
				line[0] = ' ';												
 
				//**************************************************
				// Read 3 floats from the line:  X Y Z and store them in the corresponding place in _vertices
				// In order to read the floats in one shot from string you can use sscanf
				//**************************************************
				point3d p;
				sscanf(line.c_str(),"%f %f %f ",
					&p.x,
					&p.y,
					&p.z);

				_v.push_back( p );
				_nv.push_back(vec3d());
				
				//**************************************************
				// This is for the 2nd part of the exercise: update the bounding box
				// For the very first vertex read, initialize the bb accordingly
				//**************************************************
				if( _v.size() == 1)
				{
					//**************************************************
					// Case of the very first vertex read
					//**************************************************
					_bb.set(p);
				}
				else
				{
					//**************************************************
					// otherwise check the vertex against the bounding box and in case update it
					//**************************************************
					_bb.add(p);
				}
					
			}
			//**************************************************
			// If the first character is a 'f'...
			//**************************************************
			if (line.c_str()[0] == 'f')										
			{

				triangleIndex t;
				assert(parseFaceString(line, t));

//				// Set first character to ' '. This will allow us to use sscanf
//				line[0] = ' ';
 
//				// this contains temporary the indices of the vertices
//				triangleIndex t;
				
//				//**************************************************
//				// Read 3 integers from the line:  idx1 idx2 idx3 and store them in the corresponding place in vertexIdx
//				// In order to read the 3 integers in one shot from string you can use sscanf
//				//**************************************************
//				sscanf(line.c_str(),"%hu %hu %hu",								// Read integers from the line:  f 1 2 3
//					&t.v1,										// First point of our triangle. This is an
//					&t.v2,										// pointer to our vertexBuffer list
//					&t.v3 );										// each point represents an X,Y,Z.


//				PRINTVAR(t2);
//				PRINTVAR(t);

//				assert(t2==t);

//				t=t2;

				//**************************************************
				// correct the indices: OBJ starts counting from 1, in C the arrays starts at 0...
				//**************************************************
				t -= 1;

				 //cout << vertexNumber[0] << " " << vertexNumber[1] << " " << vertexNumber[2] << " "  << endl;
				_indices.push_back( t );
				
								
				//*********************************************************************
				//  Calculate the normal of the triangles, it will be the same for each vertex
				//*********************************************************************
				vec3d norm;
				
				//*********************************************************************
				//  compute the normal for the 3 vertices we just added
				//*********************************************************************
				computeNormal( _v[ t.v1], _v[t.v2], _v[t.v3], norm );

//				PRINTVAR(angleAtVertex(_v[ t.v1], _v[t.v2], _v[t.v3] ));
//				PRINTVAR(angleAtVertex(_v[ t.v2], _v[t.v1], _v[t.v3] ));
//				PRINTVAR(angleAtVertex(_v[ t.v3], _v[t.v1], _v[t.v2] ));
//				_nv[t.v1] += (vec3d(norm) * angleAtVertex(_v[ t.v1], _v[t.v2], _v[t.v3]));
//				_nv[t.v2] += (vec3d(norm) * angleAtVertex(_v[ t.v2], _v[t.v1], _v[t.v3]));
//				_nv[t.v3] += (vec3d(norm) * angleAtVertex(_v[ t.v3], _v[t.v1], _v[t.v2]));
				_nv[t.v1] += ( vec3d(norm) );
				_nv[t.v2] += ( vec3d(norm) );
				_nv[t.v3] += ( vec3d(norm) );
				
			}
		}

		cout << "Found :\n\tNumber of triangles (_indices) "<< _indices.size() << "\n\tNumber of Vertices: " << _v.size() << "\n\tNumber of Normals: " << _nv.size() << endl;
		PRINTVAR( _indices );
		PRINTVAR( _v );
		PRINTVAR( _nv );
		
		// normalize the sum of normals
		for(int i=0; i<_nv.size(); _nv[i++].normalize());
		PRINTVAR( _nv );
		
		// Close OBJ file
		objFile.close();

	}
	else 
	{
		cout << "Unable to open file";								
	}
	

	cout << "Object loaded with "<< _v.size() << " vertices and "  << _indices.size() << " faces" << endl;
	cout << "Bounding box : pmax=" << _bb.pmax << "  pmin=" << _bb.pmin << endl;
	return 0;
}


/**
 * Computes the angle at vertex baseV formed by the edges connecting it with the
 * vertices v1 and v2 respectively, ie the baseV-v1 and baseV-v2 edges
 * @brief Computes the angle at vertex
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
//	PRINTVAR((e1).dot(e2) / (e1.norm()*e2.norm()));
//	PRINTVAR( acos( (e1).dot(e2) / (e1.norm()*e2.norm()) ));
	//safe acos...
	if(fabs( (e1).dot(e2) / (e1.norm()*e2.norm()) ) >= 1.0f )
	{
		 cerr << "warning: using safe acos" << endl;
		 return (acos(1.0f));
	}
	else
	{
		return ( acos( (e1).dot(e2) / (e1.norm()*e2.norm()) ) );
	}
}

/**
 * @brief ObjModel::subdivision
 */
void ObjModel::linearSubdivision()
{
	// copy the data in new arrays
	_subVert = _v;
	_subIdx = _indices;
	_subNorm = _nv;

	PRINTVAR(_subVert);
	PRINTVAR(_v);
	// create a list of the new vertices creates witht the reference to the edge
	EdgeList newVertices;

    // for each triangle
	for(int i = 0; i < _indices.size(); ++i )
	{
		// get the indices of the triangle vertices
		idxtype v1 = _indices[i].v1;
		idxtype v2 = _indices[i].v2;
		idxtype v3 = _indices[i].v3;

		// for each edge get the index of the vertex of the midpoint
		// if it does not exist it will be generated.
		idxtype a = getNewVertex( edge(v1,v2), _subVert, _subNorm, newVertices );
		idxtype b = getNewVertex( edge(v2,v3), _subVert, _subNorm, newVertices );
		idxtype c = getNewVertex( edge(v3,v1), _subVert, _subNorm ,newVertices );

		// create the four new triangle
		// BE CAREFULL WITH THE VERTEX ORDER!!
		//		       v2
		//		      / \
		//		     /   \
		//		   a ----- b
		//		  /  \    / \
		//		 /     \ /    \
		//		v1 ---- c ---- v3
		//
		// the original triangle was v1-v2-v3, use the same clock-wise order for the other

		_subIdx.push_back( triangleIndex( v1, a, c ) );
		_subIdx.push_back( triangleIndex( a, b, c ) );
		_subIdx.push_back( triangleIndex( a, v2, b ) );
		_subIdx.push_back( triangleIndex( c, b, v3 ) );
	}

	PRINTVAR(newVertices);
}

/**
 * For a given edge it returns the index of the new vertex created on its middle point. If such vertex already exists it just returns the
 * its index; if it does not exist it creates it in vertList along it's normal and return the index
 * @brief ObjModel::getNewVertex
 * @param e the edge
 * @param vertList the list of vertices
 * @param normList the list of normals associated to the vertices
 * @param newVertList The list of the new vertices added so far
 * @return the index of the new vertex
 * @see EdgeList
 */
idxtype ObjModel::getNewVertex( const edge &e, vector<point3d> &vertList, vector<vec3d> &normList, EdgeList &newVertList ) const
{
//	PRINTVAR(e);
//	PRINTVAR(newVertList);
    // if the egde is not contained in the new vertex list
	if( !newVertList.contains( e ) )
    {
        // generate new index (vertex.size)
        idxtype idxnew = vertList.size();
//		PRINTVAR(vertList);
//		PRINTVAR(idxnew);
        // add the edge and index to the map
		newVertList.add( e, idxnew );
//		PRINTVAR(newVertList);
        // generate new vertex
        point3d nvert = (vertList[e.first] + vertList[e.second])*0.5;
		// append it to the list of vertices
//	cout << "new vertex created " << endl;
//	PRINTVAR(nvert);
//	PRINTVAR(e);
        vertList.push_back( nvert );
		// simple version, compute the normal as the linear combination of the
		// normals of the 2 vertices
		vec3d nnorm = (normList[e.first] + normList[e.second]);
		nnorm.normalize();
		normList.push_back( nnorm );

		return idxnew;

    }
    // else
    {
        // get the index of the vertex
		return ( newVertList.getIndex( e ) );
    }
}


void ObjModel::render( const RenderingParameters &params ) 
{

	if(!params.subdivision )
	{
		draw( _v, _indices, _nv, params); 
	}
	else
	{
		if( _subIdx.empty() || _subNorm.empty() || _subVert.empty() )
		{
			linearSubdivision();
		}
		draw( _subVert, _subIdx, _subNorm, params);
	}
}
void ObjModel::draw( const vector<point3d> &vertices, const vector<triangleIndex> &indices, vector<point3d> &vertexNormals, const RenderingParameters &params ) const
{
	if( params.solid )
	{
		drawSolid( vertices, indices, vertexNormals, params );
	}
	if( params.wireframe )
	{
		drawWireframe( vertices, indices, params );
	}
	
}
void ObjModel::drawSolid( const vector<point3d> &vertices, const vector<triangleIndex> &indices, vector<point3d> &vertexNormals, const RenderingParameters &params ) const
{
	 if (params.useIndexRendering)
	 {
		 indexDraw( vertices, indices, vertexNormals, params );
	 }
	 else
	 {
		 flatDraw( vertices, indices, params );
	 }
}

void ObjModel::drawWireframe( const vector<point3d> &vertices, const vector<triangleIndex> &indices, const RenderingParameters &params ) const
{
	glDisable(GL_LIGHTING);
	if( params.solid )
	{
		glColor3f(0,0,0);
		glLineWidth(2);
	}
	else
	{
		glColor3f(.8,.8,.8);
		glLineWidth(.21);
	}
	
	for(int i=0; i<indices.size(); i++)
	{
		glBegin(GL_LINE_LOOP);
			glVertex3fv((float*)&vertices[indices[i].v1]);

			glVertex3fv((float*)&vertices[indices[i].v2]);

			glVertex3fv((float*)&vertices[indices[i].v3]);
		glEnd();
	}
	glEnable(GL_LIGHTING);
}

void ObjModel::flatDraw( const vector<point3d> &vertices, const vector<triangleIndex> &indices, const RenderingParameters &params ) const
{
	if(params.smooth)
	{
		glShadeModel( GL_SMOOTH );
	}
	else
	{
		glShadeModel( GL_FLAT );
	}

	// for each triangle draw the vertices and the normals
	for(int i=0; i<indices.size(); i++)
	{
		glBegin(GL_TRIANGLES);
			//compute the normal of the triangle
			vec3d n;
			computeNormal( vertices[indices[i].v1], vertices[indices[i].v2], vertices[indices[i].v3], n);
			glNormal3fv((float*)&n);

			glVertex3fv((float*)&vertices[indices[i].v1]);

			glVertex3fv((float*)&vertices[indices[i].v2]);

			glVertex3fv((float*)&vertices[indices[i].v3]);

		glEnd();
	}
}

void ObjModel::indexDraw( const vector<point3d> &vertices, const vector<triangleIndex> &indices, vector<point3d> &vertexNormals, const RenderingParameters &params ) const
{
	if(params.smooth)
	{
		glShadeModel( GL_SMOOTH );
	}
	else
	{
		glShadeModel( GL_FLAT );
	}
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
	glEnableClientState(GL_VERTEX_ARRAY);

	//****************************************
	// Normal pointer to normal array
	//****************************************
	glNormalPointer(GL_FLOAT, 0, (float*)&vertexNormals[0]);

	//****************************************
	// Index pointer to normal array
	//****************************************
	glVertexPointer(COORD_PER_VERTEX, GL_FLOAT, 0, (float*)&vertices[0]);

	//****************************************
	// Draw the triangles
	//****************************************
	glDrawElements(GL_TRIANGLES, indices.size()*VERTICES_PER_TRIANGLE, GL_UNSIGNED_INT, (idxtype*)&indices[0]);

	//****************************************
	// Disable vertex arrays
	//****************************************
	glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays

	//****************************************
	// Disable normal arrays
	//****************************************
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
	if( !_v.empty() && !_indices.empty() )
	{     
		//****************************************
		// calculate model width, height, and 
		// depth using the bounding box
		//****************************************
		float w,h,d;
		w = fabs(_bb.pmax.x - _bb.pmin.x);
		h = fabs(_bb.pmax.y - _bb.pmin.y);
		d = fabs(_bb.pmax.z - _bb.pmin.z);
		
cout << "size: w: " << w << " h " << h << " d " << d <<  endl;
		//****************************************
		// calculate center of the bounding box of the model 
		//****************************************
		point3d c = (_bb.pmax + _bb.pmin) * 0.5;
		
		//****************************************
		// calculate the unitizing scale factor as the 
		// maximum of the 3 dimensions
		//****************************************
		float scale = 2.0 / std::max(std::max(w, h), d);
		
cout << "scale: " << scale << " cx " << c.x << " cy " << c.y << " cz " << c.z << endl;

		// translate each vertex wrt to the center and then apply the scaling to the coordinate
		for (int i = 0; i < _v.size(); i++)
		{
			//****************************************
			// translate the vertex
			//****************************************
			_v[i].translate(-c.x, -c.y, -c.z);

			//****************************************
			// apply the scaling
			//****************************************
			_v[i].scale(scale);
			
		}


		//****************************************
		// update the bounding box, ie translate and scale the 6 coordinates
		//****************************************
		_bb.pmax = (_bb.pmax - c) * scale;
		_bb.pmin = (_bb.pmin - c) * scale;

		
		cout << "New bounding box : pmax=" << _bb.pmax << "  pmin=" << _bb.pmin << endl;
		
		return scale;
		
	}
	
	return 0;
}

void ObjModel::release( )
{
}

/**
 * @brief ObjModel::parseFaceString
 * @param toParse
 * @param out
 * @return
 */
bool ObjModel::parseFaceString( const string &toParse, triangleIndex &out) const
{
	if (toParse.c_str()[0] == 'f')
	{
		idxtype a;
		// now check the different formats: %d, %d//%d, %d/%d, %d/%d/%d
		if (strstr(toParse.c_str(), "//"))
		{
			// v//n
			return ( sscanf(toParse.c_str(), "f %u//%u %u//%u %u//%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a) == 6 );
		}
		else if (sscanf(toParse.c_str(), "f %u/%u/%u", &a, &a, &a) == 3)
		{
			// v/t/n
			return ( sscanf(toParse.c_str(), "f %u/%u/%u %u/%u/%u %u/%u/%u", &(out.v1), &a, &a, &(out.v2), &a, &a, &(out.v3), &a, &a) == 9 );
		}
		else if (sscanf(toParse.c_str(), "f %u/%u", &a, &a) == 2)
		{
			// v/t .
			return ( sscanf(toParse.c_str(), "f %u/%u %u/%u %u/%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a) == 6 );
		}
		else
		{
			// v
			sscanf(toParse.c_str(), "f %u %u %u", &(out.v1), &(out.v2), &(out.v3));
			PRINTVAR(out);
			return ( sscanf(toParse.c_str(), "f %u %u %u", &(out.v1), &(out.v2), &(out.v3)) == 3 );
		}
	}
	else
	{
		return false;
	}
}





//*****************************************************************************
//*						DEPRECATED FUNCTIONS
//*****************************************************************************

// to be deprecated
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

// to be deprecated
void ObjModel::drawWireframe() const
{

	drawWireframe(_v, _indices, RenderingParameters() );
	
}

// to be deprecated
void ObjModel::indexDraw() const
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
	glEnableClientState(GL_VERTEX_ARRAY);

	//****************************************
	// Normal pointer to normal array
	//****************************************
	glNormalPointer(GL_FLOAT, 0, (float*)&_nv[0]);

	//****************************************
	// Index pointer to normal array
	//****************************************
	glVertexPointer(COORD_PER_VERTEX, GL_FLOAT, 0, (float*)&_v[0]);

	//****************************************
	// Draw the triangles
	//****************************************
	glDrawElements(GL_TRIANGLES, _indices.size()*VERTICES_PER_TRIANGLE, GL_UNSIGNED_INT, (idxtype*)&_indices[0]);

	//****************************************
	// Disable vertex arrays
	//****************************************
	glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays

	//****************************************
	// Disable normal arrays
	//****************************************
	glDisableClientState(GL_NORMAL_ARRAY);
}

// to be deprecated
void ObjModel::drawSubdivision()
{
	if( _subIdx.empty() || _subNorm.empty() || _subVert.empty() )
	{
		linearSubdivision();
	}

	glShadeModel( GL_SMOOTH );

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);

	glNormalPointer(GL_FLOAT, 0, (float*)&_subNorm[0]);
	glVertexPointer(COORD_PER_VERTEX, GL_FLOAT, 0, (float*)&_subVert[0]);

	glDrawElements(GL_TRIANGLES, _subIdx.size()*VERTICES_PER_TRIANGLE, GL_UNSIGNED_SHORT, (idxtype*)&_subIdx[0]);


	glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);

	drawWireframe(_subVert, _subIdx, RenderingParameters() );

}