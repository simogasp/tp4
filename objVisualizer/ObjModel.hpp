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



/**
 */
typedef struct BoundingBox
{
	point3d pmax;
	point3d pmin;

	BoundingBox() :
	pmax(point3d()),
	pmin(point3d()) { }

	BoundingBox(const point3d &p) :
	pmax(p),
	pmin(p) { }

	/**
	 * Add a point to the bounding box. Its coordinates are taken into account 
	 * and the limits of the bounding box updated accordingly
	 * @param p the point to add
	 */
	void add(const point3d &p)
	{
		pmax.max(p);
		pmin.min(p);
	}

	/**
	 * Set the bounding box to the given point
	 * @param p the point
	 */
	void set(const point3d &p)
	{
		pmax = p;
		pmin = p;
	}

} BoundingBox;

enum obj2render
{
	ORIGINAL, SUBDIVISION
};

typedef struct RenderingParameters
{
	bool wireframe; //!< wireframe on/off
	bool solid; //!< draw the mesh on/off
	bool useIndexRendering; //!< use opengl drawElements on/off
	bool subdivision; //!< subdivision on/off
	bool smooth; //!< GL_SMOOTH on/off
	bool normals; //!< show normals on/off

	RenderingParameters() :
	wireframe(true),
	solid(true),
	useIndexRendering(false),
	smooth(false),
	subdivision(false),
	normals(true) { }

} RenderingParameters;

/**
 * The class containing and managing the 3D model 
 */
class ObjModel
{
private:

	std::vector<triangleIndex> _indices; //!< Stores the vertex indices for the triangles
	std::vector<point3d> _v; //!< Stores the vertices
	std::vector<vec3d> _nv; //!< Stores the normals for the triangles

	// Subdivision
	std::vector<triangleIndex> _subIdx; //!< Stores the vertex indices for the triangles
	std::vector<point3d> _subVert; //!< Stores the vertices
	std::vector<vec3d> _subNorm; //!< Stores the normals for the triangles

	BoundingBox _bb; //!< the current bounding box of the model

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
	void computeNormal(const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm) const;

	/**
	 * Computes the angle at vertex baseV formed by the edges connecting it with the
	 * vertices v1 and v2 respectively, ie the baseV-v1 and baseV-v2 edges
	 * @brief Computes the angle at vertex
	 * @param baseV the vertex at which to compute the angle
	 * @param v1 the other vertex of the first edge baseV-v1
	 * @param v2 the other vertex of the second edge baseV-v2
	 * @return the angle in radiants
	 */
	float angleAtVertex(const point3d& v1, const point3d& v2, const point3d& v3) const;

	/**
	 *  Loads the model from file
	 * @param filename the OBJ file
	 * @return 
	 */
	int load(char *filename);

	/**
	 * Render the model according to the provided parameters
	 * @param params The rendering parameters
	 */
	void render(const RenderingParameters &params = RenderingParameters());

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
	void draw(const std::vector<point3d> &vertices, const std::vector<triangleIndex> &indices, std::vector<vec3d> &vertexNormals, const RenderingParameters &params) const;
	void drawSolid(const std::vector<point3d> &vertices, const std::vector<triangleIndex> &indices, std::vector<vec3d> &vertexNormals, const RenderingParameters &params) const;
	void drawWireframe(const std::vector<point3d> &vertices, const std::vector<triangleIndex> &indices, const RenderingParameters &params) const;
	void indexDraw(const std::vector<point3d> &vertices, const std::vector<triangleIndex> &indices, std::vector<vec3d> &vertexNormals, const RenderingParameters &params) const;
	void flatDraw(const std::vector<point3d> &vertices, const std::vector<triangleIndex> &indices, const RenderingParameters &params) const;
	void drawNormals(const std::vector<point3d> &vertices, std::vector<vec3d> &vertexNormals) const;
	/////////////////////////////
	// DEPRECATED METHODS
	DEPRECATED(void drawSubdivision());
	DEPRECATED(void indexDraw() const);
	DEPRECATED(void flatDraw() const);
	DEPRECATED(void drawWireframe() const);



private:

	/**
	 * Compute the subdivision of the input mesh by applying one step of the Loop algorithm
	 * 
	 * @param[in] origVert The list of the input vertices
	 * @param[in] origMesh The input mesh (the vertex indices for each face/triangle)
	 * @param[out] destVert The list of the new vertices for the subdivided mesh
	 * @param[out] destMesh The new subdivided mesh (the vertex indices for each face/triangle)
	 * @param[out] destNorm The new list of normals for each new vertex of the subdivided mesh
	 */
	void loopSubdivision(const std::vector<point3d> &origVert, const std::vector<triangleIndex> &origMesh, std::vector<point3d> &destVert, std::vector<triangleIndex> &destMesh, std::vector<vec3d> &destNorm) const;


	/**
	 * For a given edge it returns the index of the new vertex created on its middle point. If such vertex already exists it just returns the
	 * its index; if it does not exist it creates it in vertList along it's normal and return the index
	 * @brief ObjModel::getNewVertex
	 * @param e the edge
	 * @param currFace the current triangle containing the edge e
	 * @param vertList the list of vertices
	 * @param indices the list of triangles
	 * @param normList the list of normals associated to the vertices
	 * @param newVertList The list of the new vertices added so far
	 * @return the index of the new vertex
	 * @see EdgeList
	 */
	idxtype getNewVertex(const edge &e, std::vector<point3d> &vertList, const std::vector<triangleIndex> &indices, EdgeList &newVertList) const;

	/**
	 * @brief ObjModel::parseFaceString
	 * @param toParse
	 * @param out
	 * @return
	 */
	bool parseFaceString(const std::string &toParse, triangleIndex &out) const;

	DEPRECATED(void applyLoop(const triangleIndex &t, const std::vector<point3d> &orig, std::vector<size_t> &valence, std::vector<point3d> &dest) const);

};

#endif

