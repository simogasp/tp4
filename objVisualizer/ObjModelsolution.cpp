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

ObjModel::ObjModel( ) { }

/**
 * Calculate the normal of a triangular face defined by three points
 *
 * @param[in] v1 the first vertex
 * @param[in] v2 the second vertex
 * @param[in] cv3 the third vertex
 * @param[out] norm the normal
 */
void ObjModel::computeNormal( const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm ) const
{
	norm = (v1 - v2).cross( v1 - v3 );

	norm.normalize( );
}

/**
 * @brief 
 * @param filename
 * @return 
 */
int ObjModel::load( char* filename )
{

	string line;
	ifstream objFile( filename );

	// If obj file is open, continue
	if ( objFile.is_open( ) )
	{
		// Start reading file data
		while( !objFile.eof( ) )
		{
			//**************************************************
			// Get a line from file
			//**************************************************
			getline( objFile, line );

			//**************************************************
			// If the first character is a 'v'...
			//**************************************************
			PRINTVAR( line );
			if ( (line.c_str( )[0] == 'v') && (line.c_str( )[1] == ' ') ) // to drop all the vn and vn lines
			{
				PRINTVAR( line );
				//**************************************************
				// Read 3 floats from the line:  X Y Z and store them in the corresponding place in _vertices
				// In order to read the floats in one shot from string you can use sscanf
				//**************************************************
				point3d p;
				sscanf( line.c_str( ), "v %f %f %f ",
					 &p.x,
					 &p.y,
					 &p.z );

				_v.push_back( p );
				_nv.push_back( vec3d( ) );

				//**************************************************
				// This is for the 2nd part of the exercise: update the bounding box
				// For the very first vertex read, initialize the bb accordingly
				//**************************************************
				if ( _v.size( ) == 1 )
				{
					//**************************************************
					// Case of the very first vertex read
					//**************************************************
					_bb.set( p );
				}
				else
				{
					//**************************************************
					// otherwise check the vertex against the bounding box and in case update it
					//**************************************************
					_bb.add( p );
				}

			}
			//**************************************************
			// If the first character is a 'f'...
			//**************************************************
			if ( line.c_str( )[0] == 'f' )
			{

				triangleIndex t;
				assert( parseFaceString( line, t ) );

				//**************************************************
				// correct the indices: OBJ starts counting from 1, in C the arrays starts at 0...
				//**************************************************
				t -= 1;

				_indices.push_back( t );


				//*********************************************************************
				//  Calculate the normal of the triangles, it will be the same for each vertex
				//*********************************************************************
				vec3d norm;

				//*********************************************************************
				//  compute the normal for the 3 vertices we just added
				//*********************************************************************
				computeNormal( _v[ t.v1], _v[t.v2], _v[t.v3], norm );

				_nv[t.v1] += (vec3d( norm ) * angleAtVertex( _v[ t.v1], _v[t.v2], _v[t.v3] ));
				_nv[t.v2] += (vec3d( norm ) * angleAtVertex( _v[ t.v2], _v[t.v1], _v[t.v3] ));
				_nv[t.v3] += (vec3d( norm ) * angleAtVertex( _v[ t.v3], _v[t.v1], _v[t.v2] ));
				//				_nv[t.v1] += ( vec3d(norm) );
				//				_nv[t.v2] += ( vec3d(norm) );
				//				_nv[t.v3] += ( vec3d(norm) );

			}
		}

		cerr << "Found :\n\tNumber of triangles (_indices) " << _indices.size( ) << "\n\tNumber of Vertices: " << _v.size( ) << "\n\tNumber of Normals: " << _nv.size( ) << endl;
		PRINTVAR( _indices );
		PRINTVAR( _v );
		PRINTVAR( _nv );

		// normalize the sum of normals
		for ( int i = 0; i < _nv.size( ); _nv[i++].normalize( ) );
		PRINTVAR( _nv );

		// Close OBJ file
		objFile.close( );

	}
	else
	{
		cout << "Unable to open file";
	}


	cout << "Object loaded with " << _v.size( ) << " vertices and " << _indices.size( ) << " faces" << endl;
	cout << "Bounding box : pmax=" << _bb.pmax << "  pmin=" << _bb.pmin << endl;
	return 0;
}

/**
 * Computes the angle at vertex baseV formed by the edges connecting it with the
 * vertices v1 and v2 respectively, ie the baseV-v1 and baseV-v2 edges
 * 
 * @brief Computes the angle at vertex
 * @param baseV the vertex at which to compute the angle
 * @param v1 the other vertex of the first edge baseV-v1
 * @param v2 the other vertex of the second edge baseV-v2
 * @return the angle in radiants
 */
float ObjModel::angleAtVertex( const point3d& baseV, const point3d& v1, const point3d& v2 ) const
{
	vec3d e1 = baseV - v1;
	vec3d e2 = baseV - v2;
	//safe acos...
	if ( fabs( (e1).dot( e2 ) / (e1.norm( ) * e2.norm( )) ) >= 1.0f )
	{
		cerr << "warning: using safe acos" << endl;
		return (acos( 1.0f ));
	}
	else
	{
		return ( acos( (e1).dot( e2 ) / (e1.norm( ) * e2.norm( )) ));
	}
}

/**
 * Compute the subdivision of the input mesh by applying one step of the Loop algorithm
 * 
 * @param[in] origVert The list of the input vertices
 * @param[in] origMesh The input mesh (the vertex indices for each face/triangle)
 * @param[out] destVert The list of the new vertices for the subdivided mesh
 * @param[out] destMesh The new subdivided mesh (the vertex indices for each face/triangle)
 * @param[out] destNorm The new list of normals for each new vertex of the subdivided mesh
 */
void ObjModel::loopSubdivision( const std::vector<point3d> &origVert,			//!< the original vertices
								const std::vector<triangleIndex> &origMesh,		//!< the original mesh
								std::vector<point3d> &destVert,					//!< the new vertices
								std::vector<triangleIndex> &destMesh,			//!< the new mesh
								std::vector<vec3d> &destNorm ) const			//!< the new normals
{
	// copy the data in new arrays
	destVert = origVert;
	destMesh.clear( );

	//	PRINTVAR(destVert);
	//	PRINTVAR(origVert);
	// create a list of the new vertices creates with the reference to the edge
	EdgeList newVertices;

	// for each triangle
	for ( int i = 0; i < origMesh.size( ); ++i )
	{
		// get the indices of the triangle vertices
		idxtype v1 = origMesh[i].v1;
		idxtype v2 = origMesh[i].v2;
		idxtype v3 = origMesh[i].v3;

		// for each edge get the index of the vertex of the midpoint
		// if it does not exist it will be generated.
		idxtype a = getNewVertex( edge( v1, v2 ), destVert, origMesh, newVertices );
		idxtype b = getNewVertex( edge( v2, v3 ), destVert, origMesh, newVertices );
		idxtype c = getNewVertex( edge( v3, v1 ), destVert, origMesh, newVertices );

		// create the four new triangle
		// BE CAREFULL WITH THE VERTEX ORDER!!
		//		       v2
		//			   /\
		//		      /  \
		//		     /    \
		//		    a ---- b
		//         / \     /\
		//		  /   \   /  \
		//		 /     \ /    \
		//		v1 ---- c ---- v3
		//
		// the original triangle was v1-v2-v3, use the same clock-wise order for the other
		// hence v1-a-c, a-b-c and so on

		destMesh.push_back( triangleIndex( v1, a, c ) );
		destMesh.push_back( triangleIndex( a, b, c ) );
		destMesh.push_back( triangleIndex( a, v2, b ) );
		destMesh.push_back( triangleIndex( c, b, v3 ) );
	}

	// if initialize at 0 need to update the count before applying loop
	vector<size_t> valence( destVert.size( ), 0 );

	vector<point3d> tmp( destVert.size( ) ); // a copy

	for ( int i = 0; i < origMesh.size( ); ++i )
	{
		triangleIndex t = origMesh[i];
		
		valence[t.v1]++;
		tmp[t.v1] += (0.625f * origVert[t.v1] + 0.1875f * origVert[t.v2] + 0.1875f * origVert[t.v3]);

		valence[t.v2]++;
		tmp[t.v2] += (0.625f * origVert[t.v2] + 0.1875f * origVert[t.v1] + 0.1875f * origVert[t.v3]);

		valence[t.v3]++;
		tmp[t.v3] += (0.625f * origVert[t.v3] + 0.1875f * origVert[t.v2] + 0.1875f * origVert[t.v1]);
	}

	for ( int i = 0; i < origVert.size( ); ++i )
	{
		assert( valence[i] != 0 );
		destVert[i] = tmp[i] / valence[i];
	}
	//PRINTVAR(destVert);
	destNorm.clear( );
	destNorm = vector<vec3d>(destVert.size( ));

	//for each face
	for ( int i = 0; i < destMesh.size( ); ++i )
	{
		//*********************************************************************
		//  Calculate the normal of the triangles, it will be the same for each vertex
		//*********************************************************************
		vec3d norm;

		//*********************************************************************
		//  compute the normal for the 3 vertices we just added
		//*********************************************************************
		computeNormal( destVert[ destMesh[i].v1], destVert[destMesh[i].v2], destVert[destMesh[i].v3], norm );

		destNorm[destMesh[i].v1] += (vec3d( norm ));
		destNorm[destMesh[i].v2] += (vec3d( norm ));
		destNorm[destMesh[i].v3] += (vec3d( norm ));

	}
	for ( int i = 0; i < destNorm.size( ); destNorm[i++].normalize( ) );
	//PRINTVAR(newVertices);
}



/**
 * For a given edge it returns the index of the new vertex created on its middle point. If such vertex already exists it just returns the
 * its index; if it does not exist it creates it in vertList along it's normal and return the index
 * 
 * @param e the edge
 * @param currFace the current triangle containing the edge e
 * @param vertList the list of vertices
 * @param indices the list of triangles
 * @param normList the list of normals associated to the vertices
 * @param newVertList The list of the new vertices added so far
 * @return the index of the new vertex
 * @see EdgeList
 */
idxtype ObjModel::getNewVertex( const edge &e,
								std::vector<point3d> &vertList,
								const std::vector<triangleIndex> &indices,
								EdgeList &newVertList ) const
{
	//	PRINTVAR(e);
	//	PRINTVAR(newVertList);
	// if the egde is not contained in the new vertex list
	if ( !newVertList.contains( e ) )
	{
		// generate new index (vertex.size)
		idxtype idxnew = vertList.size( );
		// add the edge and index to the map
		newVertList.add( e, idxnew );
		// generate new vertex
		point3d nvert;
		idxtype oppV1;
		idxtype oppV2;
		// check if it is a boundary edge, ie check if there is another triangle 
		// sharing this edge and if so get the index of its "opposite" vertex
		// if it is not a boundary
		if ( !isBoundaryEdge( e, indices, oppV1, oppV2 ) )
		{
			// the new vertex is the linear combination of the two extrema of 
			//the edge V1 and V2 and the two opposite vertices oppV1 and oppV2
			// Using the loop coefficient the new vertex is
			// nvert = 3/8 (V1+V2) + 1/8(oppV1 + oppV2)
			nvert = 0.375f * (vertList[e.first] + vertList[e.second]) + 0.125f * (vertList[oppV1] + vertList[oppV2]);
		}
		else
		{
			// otherwise it is a boundary edge then the vertex is the linear combination of the 
			// two extrema
			nvert = (vertList[e.first] + vertList[e.second])*0.5;
		}
		// append it to the list of vertices
		vertList.push_back( nvert );
		return idxnew;

	}
	// else
	{
		// get the index of the vertex
		return ( newVertList.getIndex( e ));
	}
}



void ObjModel::drawWireframe( const vector<point3d> &vertices, const vector<triangleIndex> &indices, const RenderingParameters &params ) const
{
	glDisable( GL_LIGHTING );
	if ( params.solid )
	{
		glColor3f( 0, 0, 0 );
		glLineWidth( 2 );
	}
	else
	{
		glColor3f( .8, .8, .8 );
		glLineWidth( .21 );
	}

	for ( int i = 0; i < indices.size( ); i++ )
	{
		glBegin( GL_LINE_LOOP );
		glVertex3fv( (float*) &vertices[indices[i].v1] );

		glVertex3fv( (float*) &vertices[indices[i].v2] );

		glVertex3fv( (float*) &vertices[indices[i].v3] );
		glEnd( );
	}
	glEnable( GL_LIGHTING );
}

void ObjModel::flatDraw( const vector<point3d> &vertices, const vector<triangleIndex> &indices, const RenderingParameters &params ) const
{
	if ( params.smooth )
	{
		glShadeModel( GL_SMOOTH );
	}
	else
	{
		glShadeModel( GL_FLAT );
	}

	// for each triangle draw the vertices and the normals
	for ( int i = 0; i < indices.size( ); i++ )
	{
		glBegin( GL_TRIANGLES );
		//compute the normal of the triangle
		vec3d n;
		computeNormal( vertices[indices[i].v1], vertices[indices[i].v2], vertices[indices[i].v3], n );
		glNormal3fv( (float*) &n );

		glVertex3fv( (float*) &vertices[indices[i].v1] );

		glVertex3fv( (float*) &vertices[indices[i].v2] );

		glVertex3fv( (float*) &vertices[indices[i].v3] );

		glEnd( );
	}
}

void ObjModel::drawNormals( const std::vector<point3d> &vertices, std::vector<vec3d> &vertexNormals ) const
{
	glDisable( GL_LIGHTING );

	glColor3f( 0.8, 0, 0 );
	glLineWidth( 2 );

	for ( int i = 0; i < vertices.size( ); i++ )
	{
		glBegin( GL_LINES );

		vec3d newP = vertices[i] + 0.1 * vertexNormals[i];
		glVertex3fv( (float*) &vertices[i] );

		glVertex3f( newP.x, newP.y, newP.z );

		glEnd( );
	}
	glEnable( GL_LIGHTING );
}

void ObjModel::indexDraw( const vector<point3d> &vertices, const vector<triangleIndex> &indices, vector<vec3d> &vertexNormals, const RenderingParameters &params ) const
{
	if ( params.smooth )
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
	glEnableClientState( GL_VERTEX_ARRAY );

	//****************************************
	// Enable normal arrays
	//****************************************
	glEnableClientState( GL_NORMAL_ARRAY );

	//****************************************
	// Vertex Pointer to triangle array
	//****************************************
	glEnableClientState( GL_VERTEX_ARRAY );

	//****************************************
	// Normal pointer to normal array
	//****************************************
	glNormalPointer( GL_FLOAT, 0, (float*) &vertexNormals[0] );

	//****************************************
	// Index pointer to normal array
	//****************************************
	glVertexPointer( COORD_PER_VERTEX, GL_FLOAT, 0, (float*) &vertices[0] );

	//****************************************
	// Draw the triangles
	//****************************************
	glDrawElements( GL_TRIANGLES, indices.size( ) * VERTICES_PER_TRIANGLE, GL_UNSIGNED_INT, (idxtype*) & indices[0] );

	//****************************************
	// Disable vertex arrays
	//****************************************
	glDisableClientState( GL_VERTEX_ARRAY ); // disable vertex arrays

	//****************************************
	// Disable normal arrays
	//****************************************
	glDisableClientState( GL_NORMAL_ARRAY );

}

/**
* Render the model according to the provided parameters
* @param params The rendering parameters
*/
void ObjModel::render( const RenderingParameters &params )
{
	// if we need to draw the original model
	if ( !params.subdivision )
	{
		// draw it
		draw( _v, _indices, _nv, params );
		// draw the normals
		if ( params.normals )
		{
			drawNormals( _v, _nv );
		}
	}
	else
	{
		PRINTVAR(params.subdivLevel);
		PRINTVAR(_currentSubdivLevel);
		// before drawing check the current level of subdivision and the required one
		if ( ( _currentSubdivLevel == 0 ) || ( _currentSubdivLevel != params.subdivLevel ) )
		{
			// if they are different apply the missing steps: either restart from the beginning
			// if the required level is less than the current one or apply the missing
			// steps starting from the current one
			vector<point3d> tmpVert;		//!< a temporary list of vertices used in the iterations
			vector<triangleIndex> tmpMesh;	//!< a temporary mesh used in the iterations
			
			if(( _currentSubdivLevel == 0 ) || ( _currentSubdivLevel > params.subdivLevel ) )
			{
				// start from the beginning
				_currentSubdivLevel = 0;
				tmpVert = _v;
				tmpMesh = _indices;
			}
			else
			{
				// start from the current level
				tmpVert = _subVert;
				tmpMesh = _subIdx;
			}
			
			// apply the proper subdivision iterations
			for( ; _currentSubdivLevel < params.subdivLevel; ++_currentSubdivLevel)
			{
				cerr << "[Loop subdivision] iteration " << _currentSubdivLevel << endl;
				loopSubdivision( tmpVert, tmpMesh, _subVert, _subIdx, _subNorm );
				// swap unless it's the last iteration
				if( _currentSubdivLevel < ( params.subdivLevel - 1) )
				{
					tmpVert = _subVert;
					tmpMesh = _subIdx;
				}
			}
		}
		
		draw( _subVert, _subIdx, _subNorm, params );
		if ( params.normals )
		{
			drawNormals( _subVert, _subNorm );
		}
	}
}

void ObjModel::draw( const vector<point3d> &vertices, const vector<triangleIndex> &indices, vector<vec3d> &vertexNormals, const RenderingParameters &params ) const
{
	if ( params.solid )
	{
		drawSolid( vertices, indices, vertexNormals, params );
	}
	if ( params.wireframe )
	{
		drawWireframe( vertices, indices, params );
	}

}

void ObjModel::drawSolid( const vector<point3d> &vertices, const vector<triangleIndex> &indices, vector<vec3d> &vertexNormals, const RenderingParameters &params ) const
{
	if ( params.useIndexRendering )
	{
		indexDraw( vertices, indices, vertexNormals, params );
	}
	else
	{
		flatDraw( vertices, indices, params );
	}
}

#include "ObjModel.cxx"