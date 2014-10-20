/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/***************************************************************************
  Simple Class to load and draw 3D objects from OBJ files
  Using triangles and normals as static object. No texture mapping. 
  OBJ files must be triangulated!!!
 ***************************************************************************/
 
#ifndef _OBJMODEL_HPP_
#define	_OBJMODEL_HPP_

#include "core.hpp"

#include <vector>
#include <ostream>
#include <math.h>

#define VERTICES_PER_TRIANGLE 3
#define COORD_PER_VERTEX 3
#define TOTAL_FLOATS_IN_TRIANGLE (VERTICES_PER_TRIANGLE*COORD_PER_VERTEX)

typedef struct v3f point3d;
typedef struct v3f vec3d;
typedef struct tindex triangleIndex;

 /**
  */
typedef struct BoundingBox
{
	point3d pmax;
	point3d pmin;
	
	BoundingBox( ) :
		pmax( point3d() ),
		pmin( point3d() ) {}

	BoundingBox( const point3d &p) :
		pmax( p ),
		pmin( p ) {}

	void add( const point3d &p )
	{
		pmax.max(p);
		pmin.min(p);
	}

	void set( const point3d &p )
	{
		pmax = p;
		pmin = p;
	}


} BoundingBox;
 




class ObjModel
{
  public: 
	ObjModel();			
		
		/**
		 * Calculate the normal of a triangular face defined by three points
		 *
		 * @param[in] v1 the first vertex
		 * @param[in] v2 the second vertex
		 * @param[in] cv3 the third vertex
		 * @param[out] norm the normal
		 */
        void computeNormal( const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm  ) const;

		/**
		 * Computes the angle at vertex baseV formed by the edges connecting it with the
		 * vertices v1 and v2 respectively, ie the baseV-v1 and baseV-v2 edges
		 * @brief Computes the angle at vertex
		 * @param baseV the vertex at which to compute the angle
		 * @param v1 the other vertex of the first edge baseV-v1
		 * @param v2 the other vertex of the second edge baseV-v2
		 * @return the angle in radiants
		 */
        float angleAtVertex( const point3d& v1, const point3d& v2, const point3d& v3 ) const;

		/**
		 * @brief ObjModel::parseFaceString
		 * @param toParse
		 * @param out
		 * @return
		 */
		bool parseFaceString( const std::string &toParse, triangleIndex &out) const;
		
		/**
		 *  Loads the model from file
		 * @param filename the OBJ file
		 * @return 
		 */
		int load(char *filename);	
		
		void drawSubdivision();

        void indexDraw() const;

        void flatDraw() const;

		void drawWireframe() const;

		void linearSubdivision();

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
		GLushort getNewVertex( const edge &e, std::vector<point3d> &vertList, std::vector<vec3d> &normList, EdgeList &newVertList ) const;
	
		/**
		 * Release the model
		 */
		void release();				 
	
		/**
		 * It scales the model to unitary size by translating it to the origin and
		 * scaling it to fit in a unit cube around the origin.
		 * 
		 * @return the scale factor used to transform the model
		 */
		float unitizeModel();

  private:
//		float* _normals;			// Stores the normals
//		float* _triangles;			// Stores the triangles
//		float* _vertices;			// Stores the points which make the object
        std::vector<triangleIndex> _indices;	// Stores the vertex indices for the triangles
        std::vector<point3d> _v;	// Stores the vertices
		std::vector<vec3d> _nt;      // Stores the normals for the triangles @todo remove ir
        std::vector<vec3d> _nv;      // Stores the normals for the triangles

		// Subdivision
		std::vector<triangleIndex> _subIdx;	// Stores the vertex indices for the triangles
		std::vector<point3d> _subVert;		// Stores the vertices
		std::vector<vec3d> _subNorm;		// Stores the normals for the triangles
					
//		long _numVertices;          // the actual number of loaded vertices
//		long _numTriangles;         // the actual number of loaded faces
		

  private:
  private:
	  
	// contains the bounding box of the model
	BoundingBox _bb;
 
};
 
#endif

