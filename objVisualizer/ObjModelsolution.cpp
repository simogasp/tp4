/**
 * @file ObjModel.cpp
 * @author  Simone Gasparini <simone.gasparini@enseeiht.fr>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * @section DESCRIPTION
 * 
 * Simple Class to load and draw 3D objects from OBJ files
 * Using triangles and normals as static object. No texture mapping. 
 * OBJ files must be triangulated!!!
 *
 */

#include "ObjModel.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cassert>
//#include <stdlib.h>

#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;


/**
 * Load the OBJ data from file
 * @param filename The name of the OBJ file
 * @return 0 if everything went well
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
            // Get a line from file
            getline( objFile, line );

            // If the first character is a simple 'v'...
//			PRINTVAR( line );
            if ( (line.c_str( )[0] == 'v') && (line.c_str( )[1] == ' ') ) // to drop all the vn and vn lines
            {
//				PRINTVAR( line );
                // Read 3 floats from the line:  X Y Z and store them in the corresponding place in _vertices
                point3d p;
                sscanf( line.c_str( ), "v %f %f %f ",  &p.x,  &p.y,  &p.z );

                //**************************************************
                // add the new point to the list of the vertices
                // and its normal to the list of normals: for the time
                // being it is a [0, 0 ,0] normal.
                //**************************************************
                _vertices.push_back( p );  //!!
                _normals.emplace_back();  //!!

                // update the bounding box, if it is the first vertex simply
                // set the bb to it
                if ( _vertices.size( ) == 1 )
                {
                    _bb.set( p );
                }
                else
                {
                    // otherwise add the point
                    _bb.add( p );
                }

            }
            // If the first character is a 'f'...
            if ( line.c_str( )[0] == 'f' )
            {

                face t;
                assert( parseFaceString( line, t ) );

                //**************************************************
                // correct the indices: OBJ starts counting from 1, in C the arrays starts at 0...
                //**************************************************
                t -= 1;  //!!

                //**************************************************
                // add it to the mesh
                //**************************************************
                _mesh.push_back( t );  //!!

                //*********************************************************************
                //  Compute the normal of the face
                //*********************************************************************
                vec3d norm;  //!!
                computeNormal( _vertices[ t.v1], _vertices[t.v2], _vertices[t.v3], norm );  //!!

                //*********************************************************************
                // Sum the normal of the face to each vertex normal
                //*********************************************************************

                _normals[t.v1] += (vec3d( norm ) * angleAtVertex( _vertices[ t.v1], _vertices[t.v2], _vertices[t.v3] ));  //!!
                _normals[t.v2] += (vec3d( norm ) * angleAtVertex( _vertices[ t.v2], _vertices[t.v1], _vertices[t.v3] ));  //!!
                _normals[t.v3] += (vec3d( norm ) * angleAtVertex( _vertices[ t.v3], _vertices[t.v1], _vertices[t.v2] ));  //!!

            }
        }

        cerr << "Found :\n\tNumber of triangles (_indices) " << _mesh.size( ) << "\n\tNumber of Vertices: " << _vertices.size( ) << "\n\tNumber of Normals: " << _normals.size( ) << endl;
//		PRINTVAR( _mesh );
//		PRINTVAR( _vertices );
//		PRINTVAR( _normals );

        //*********************************************************************
        // normalize the normals of each vertex
        //*********************************************************************
        for(auto &normal : _normals)  //<!!
        {
            normal.normalize();
        }  //>!!

//		PRINTVAR( _normals );

        // Close OBJ file
        objFile.close( );

    }
    else
    {
        cout << "Unable to open file";
    }


    cout << "Object loaded with " << _vertices.size( ) << " vertices and " << _mesh.size( ) << " faces" << endl;
    cout << "Bounding box : pmax=" << _bb.pmax << "  pmin=" << _bb.pmin << endl;
    return 0;
}

/**
 * Draw the wireframe of the model
 * 
 * @param vertices The list of vertices
 * @param mesh The mesh as a list of faces, each face is a tripleIndex of vertex indices 
 * @param params The rendering parameters
 */
void ObjModel::drawWireframe( const std::vector<point3d> &vertices, const std::vector<face> &mesh, const RenderingParameters &params ) const
{
    //**************************************************
    // we first need to disable the lighting in order to
    // draw colored segments
    //**************************************************
    glDisable( GL_LIGHTING );

    // if we are displaying the object with colored faces
    if ( params.solid )
    {
        // use black ticker lines
        glColor3f( 0, 0, 0 );
        glLineWidth( 2 );
    }
    else
    {
        // otherwise use white thinner lines for wireframe only
        glColor3f( .8, .8, .8 );
        glLineWidth( .21 );
    }

    //**************************************************
    // for each face of the mesh...
    //**************************************************
    for(const auto &face : mesh)  //!!
    {
        //**************************************************
        // draw the contour of the face as a  GL_LINE_LOOP
        //**************************************************
        glBegin( GL_LINE_LOOP );  //<!!
            glVertex3fv( (float*) &vertices[face.v1] );

            glVertex3fv( (float*) &vertices[face.v2] );

            glVertex3fv( (float*) &vertices[face.v3] );
        glEnd( );  //>!!
    }

    //**************************************************
    // re-enable the lighting
    //**************************************************
    glEnable( GL_LIGHTING );
}


/**
 * Calculate the normal of a triangular face defined by three points
 *
 * @param[in] v1 the first vertex
 * @param[in] v2 the second vertex
 * @param[in] v3 the third vertex
 * @param[out] norm the normal
 */
void ObjModel::computeNormal( const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm ) const
{
    //**************************************************
    // compute the cross product between two edges of the triangular face
    //**************************************************
    norm = (v1 - v2).cross( v1 - v3 );  //!!

    //**************************************************
    // remember to normalize the result
    //**************************************************
    norm.normalize( );  //!!
}


/**
 * Draw the faces using the computed normal of each face
 * 
 * @param vertices The list of vertices
 * @param mesh The list of face, each face containing the indices of the vertices
 * @param params The rendering parameters
 */
void ObjModel::drawFlatFaces( const std::vector<point3d> &vertices, const std::vector<face> &mesh, const RenderingParameters &params ) const
{
    // shading model to use
    if ( params.smooth )
    {
        glShadeModel( GL_SMOOTH );
    }
    else
    {
        glShadeModel( GL_FLAT );
    }

    //**************************************************
    // for each face
    //**************************************************
    for(const auto &face : mesh)  //!!
    {
        //**************************************************
        // Compute the normal to the face and then draw the
        // faces as GL_TRIANGLES assigning the proper normal
        //**************************************************
        glBegin( GL_TRIANGLES );  //<!!

            vec3d n; //the normal of the face
            computeNormal( vertices[face.v1], vertices[face.v2], vertices[face.v3], n );
            glNormal3fv( (float*) &n );

            glVertex3fv( (float*) &vertices[face.v1] );

            glVertex3fv( (float*) &vertices[face.v2] );

            glVertex3fv( (float*) &vertices[face.v3] );

        glEnd( );  //>!!
    }
}



/**
 * Draw the model using the vertex indices
 * 
 * @param vertices The vertices 
 * @param indices The list of the faces, each face containing the 3 indices of the vertices
 * @param vertexNormals The list of normals associated to each vertex
 * @param params The rendering parameters
 */
void ObjModel::drawSmoothFaces( const std::vector<point3d> &vertices, 
                                const std::vector<face> &mesh,
                                std::vector<vec3d> &vertexNormals,
                                const RenderingParameters &params ) const
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
    glEnableClientState( GL_VERTEX_ARRAY );  //!!

    //****************************************
    // Enable normal arrays
    //****************************************
    glEnableClientState( GL_NORMAL_ARRAY );  //!!

    //****************************************
    // Normal pointer to normal array
    //****************************************
    glNormalPointer( GL_FLOAT, 0, (float*) &vertexNormals[0] );  //!!

    //****************************************
    // Vertex pointer to Vertex array
    //****************************************
    glVertexPointer( COORD_PER_VERTEX, GL_FLOAT, 0, (float*) &vertices[0] );  //!!

    //****************************************
    // Draw the faces
    //****************************************
    glDrawElements( GL_TRIANGLES, mesh.size( ) * VERTICES_PER_TRIANGLE, GL_UNSIGNED_INT, (idxtype*) & mesh[0] );  //!!

    //****************************************
    // Disable vertex arrays
    //****************************************
    glDisableClientState( GL_VERTEX_ARRAY );  //!!

    //****************************************
    // Disable normal arrays
    //****************************************
    glDisableClientState( GL_NORMAL_ARRAY );  //!!

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
                                const std::vector<face> &origMesh,		        //!< the original mesh
                                std::vector<point3d> &destVert,					//!< the new vertices
                                std::vector<face> &destMesh,			        //!< the new mesh
                                std::vector<vec3d> &destNorm ) const			//!< the new normals
{
    // copy the original vertices in destVert
    destVert = origVert;

    // start fresh with the new mesh
    destMesh.clear( );

    //	PRINTVAR(destVert);
    //	PRINTVAR(origVert);

    // create a list of the new vertices creates with the reference to the edge
    EdgeList newVertices;

    //*********************************************************************
    // for each face
    //*********************************************************************
    for(const auto &f : origMesh)  //!!
    {
        //*********************************************************************
        // get the indices of the triangle vertices
        //*********************************************************************
        idxtype v1 = f.v1;  //<!!
        idxtype v2 = f.v2;
        idxtype v3 = f.v3;  //>!!

        //*********************************************************************
        // for each edge get the index of the vertex of the midpoint using getNewVertex
        //*********************************************************************
        idxtype a = getNewVertex( edge( v1, v2 ), destVert, origMesh, newVertices );  //<!!
        idxtype b = getNewVertex( edge( v2, v3 ), destVert, origMesh, newVertices );
        idxtype c = getNewVertex( edge( v3, v1 ), destVert, origMesh, newVertices );  //>!!

        //*********************************************************************
        // create the four new triangles
        // BE CAREFUL WITH THE VERTEX ORDER!!
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
        //*********************************************************************

        destMesh.emplace_back(v1, a, c);  //<!!
        destMesh.emplace_back(a, b, c);
        destMesh.emplace_back(a, v2, b);
        destMesh.emplace_back(c, b, v3);  //>!!
    }

    //*********************************************************************
    // Update each "old" vertex using the Loop coefficients. A smart way to do
    // so is to think in terms of faces than the single vertex: for each face
    // we update each of the 3 vertices using the Loop formula wrt the other 2 and
    // sum it to a temporary copy tmp of the vertices (which is initialized to
    // [0 0 0] at the beginning). We also keep a record of the occurrence of each vertex.
    // At then end, to get the final vertices we just need to divide each vertex
    // in tmp by its occurrence
    //*********************************************************************

    // A list containing the occurrence of each vertex
    vector<size_t> occurrences( origVert.size( ), 0 );

    // A list of the same size as origVert with all the elements initialized to [0 0 0]
    vector<point3d> tmp( origVert.size( ) );

    //*********************************************************************
    // for each face
    //*********************************************************************
    for(const auto &face : origMesh)  //!!
    {
        //*********************************************************************
        // consider each of the 3 vertices:
        // 1) increment its occurrence
        // 2) apply Loop update wrt the other 2 vertices of the face
        // BE CAREFUL WITH THE COEFFICIENT OF THE OTHER 2 VERTICES!... consider
        // how many times each vertex is summed in the general case...
        //*********************************************************************

        occurrences[face.v1]++;  //<!!
        tmp[face.v1] += (0.625f * origVert[face.v1] + 0.1875f * origVert[face.v2] + 0.1875f * origVert[face.v3]);

        occurrences[face.v2]++;
        tmp[face.v2] += (0.625f * origVert[face.v2] + 0.1875f * origVert[face.v1] + 0.1875f * origVert[face.v3]);

        occurrences[face.v3]++;
        tmp[face.v3] += (0.625f * origVert[face.v3] + 0.1875f * origVert[face.v2] + 0.1875f * origVert[face.v1]);  //>!!
    }

    //*********************************************************************
    //  To obtain the new vertices, divide each vertex by its occurrence value
    //*********************************************************************
    for ( size_t i = 0; i < origVert.size( ); ++i )  //!!
    {
        assert( occurrences[i] != 0 );  //??
        destVert[i] = tmp[i] / occurrences[i];  //!!
    }
    //PRINTVAR(destVert);

    // redo the normals, reset and create a list of normals of the same size as
    // the vertices, each normal set to [0 0 0]
    destNorm.clear( );
    destNorm = vector<vec3d>(destVert.size( ));

    //*********************************************************************
    //  Recompute the normals for each face
    //*********************************************************************
    for(const auto &face : destMesh)  //!!
    {
        //*********************************************************************
        //  Calculate the normal of the triangles, it will be the same for each vertex
        //*********************************************************************
        vec3d norm;  //!!
        computeNormal( destVert[face.v1], destVert[face.v2], destVert[face.v3], norm );  //!!

        //*********************************************************************
        // Sum the normal of the face to each vertex normal using the angleAtVertex as weight
        //*********************************************************************
        destNorm[face.v1] += (vec3d( norm ) * angleAtVertex( destVert[face.v1], destVert[face.v2], destVert[face.v3] ));  //<!!
        destNorm[face.v2] += (vec3d( norm ) * angleAtVertex( destVert[face.v2], destVert[face.v3], destVert[face.v1] ));
        destNorm[face.v3] += (vec3d( norm ) * angleAtVertex( destVert[face.v3], destVert[face.v1], destVert[face.v2] ));  //>!!

    }
    //*********************************************************************
    // normalize the normals of each vertex
    //*********************************************************************
    for(auto &n : destNorm)  //<!!
    {
        n.normalize( );
    }  //>!!

}

/**
 * For a given edge it returns the index of the new vertex created on its middle point. 
 * If such vertex already exists it just returns the its index; if it does not exist 
 * it creates it in vertList along it's normal and return the index
 * 
 * @param[in] e the edge
 * @param[in,out] vertList the list of vertices
 * @param[in] mesh the list of triangles
 * @param[in,out] newVertList The list of the new vertices added so far
 * @return the index of the new vertex or the one that has been already created for that edge
 * @see EdgeList
 */
idxtype ObjModel::getNewVertex( const edge &e,
                                std::vector<point3d> &vertList,
                                const std::vector<face> &mesh,
                                EdgeList &newVertList ) const
{
    //	PRINTVAR(e);
    //	PRINTVAR(newVertList);

    //*********************************************************************
    // if the egde is NOT contained in the new vertex list (see EdgeList.contains() method)
    //*********************************************************************
    if ( !newVertList.contains( e ) )  //!!
    {
        //*********************************************************************
        // generate new index (vertex.size)
        //*********************************************************************
        idxtype idxnew = static_cast<idxtype>(vertList.size( ));  //!!

        //*********************************************************************
        // add the edge and index to the newVertList
        //*********************************************************************
        newVertList.add( e, idxnew );  //!!

        // generate new vertex
        point3d nvert;		//!< this will contain the new vertex
        idxtype oppV1;		//!< the index of the first "opposite" vertex
        idxtype oppV2;		//!< the index of the second "opposite" vertex (if it exists)

        //*********************************************************************
        // check if it is a boundary edge, ie check if there is another triangle
        // sharing this edge and if so get the index of its "opposite" vertex
        //*********************************************************************
        if ( !isBoundaryEdge( e, mesh, oppV1, oppV2 ) )  //!!
        {
            // if it is not a boundary edge create the new vertex

            //*********************************************************************
            // the new vertex is the linear combination of the two extrema of
            // the edge V1 and V2 and the two opposite vertices oppV1 and oppV2
            // Using the loop coefficient the new vertex is
            // nvert = 3/8 (V1+V2) + 1/8(oppV1 + oppV2)
            //
            // REMEMBER THAT IN THE CODE OPPV1 AND OPPV2 ARE INDICES, NOT VERTICES!!!
            //*********************************************************************

            nvert = 0.375f * ( vertList[e.first] + vertList[e.second] ) + 0.125f * (vertList[oppV1] + vertList[oppV2]);  //!!
        }
        else  //??
        {
            //*********************************************************************
            // otherwise it is a boundary edge then the vertex is the linear combination of the
            // two extrema
            //*********************************************************************
            nvert = 0.5 * (vertList[e.first] + vertList[e.second]);  //!!
        }
        //*********************************************************************
        // append the new vertex to the list of vertices
        //*********************************************************************
        vertList.push_back( nvert );  //!!

        //*********************************************************************
        // return the index of the new vertex
        //*********************************************************************
        return idxnew;  //!!

    }
    else  //??
    // else we don't need to do anything, just return the associated index of the
    // already existing vertex
    {
        //*********************************************************************
        // get and return the index of the vertex
        //*********************************************************************
        return ( newVertList.getIndex( e ));  //!!
    }
}



#include "ObjModel.cxx"
