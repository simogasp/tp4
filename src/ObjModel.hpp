/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#pragma once

#include "core.hpp"
#include "objReader.hpp"
#include "rendering.hpp"

#include <cmath>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

/**
 * The class containing and managing the 3D model 
 */
class ObjModel
{
private:

    /// Stores the vertex indices for the triangles
    std::vector<face> _mesh{};
    /// Stores the vertices
    std::vector<point3d> _vertices{};
    /// Stores the normals for the triangles
    std::vector<vec3d> _normals{};

    // Subdivision
    /// Stores the vertex indices for the triangles
    std::vector<face> _subMesh{};
    /// Stores the vertices
    std::vector<point3d> _subVert{};
    /// Stores the normals for the triangles
    std::vector<vec3d> _subNorm{};

    /// the current bounding box of the model
    BoundingBox _bb{};

    /// the current subdivision level
    unsigned short _currentSubdivLevel{};   

public:
    ObjModel() = default;

    /**
     * Load the OBJ data from file
      * @param[in] filename The name of the OBJ file
      * @return true if everything went well, false otherwise
     */
    bool load(const std::string& filename);

    /**
     * Render the model according to the provided parameters
     * @param params The rendering parameters
     */
    void render(const RenderingParameters &params = RenderingParameters());


    /**
     * It scales the model to unitary size by translating it to the origin and
     * scaling it to fit in a unit cube around the origin.
     *
     * @return the scale factor used to transform the model
     */
    float unitizeModel();


private:
    /**
    * Draw the model
    *
    * @param[in] vertices list of vertices
    * @param[in] indices list of faces
    * @param[in] vertexNormals list of normals
    * @param[in] params Rendering parameters
    */
    void draw(const std::vector<point3d> &vertices, const std::vector<face> &indices, std::vector<vec3d> &vertexNormals, const RenderingParameters &params) const;


    /////////////////////////////
    // DEPRECATED METHODS
    [[deprecated]] void drawSubdivision();
    [[deprecated]] void indexDraw() const;
    [[deprecated]] void flatDraw() const;
    [[deprecated]] void drawWireframe() const;



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
    void loopSubdivision(const std::vector<point3d> &origVert, const std::vector<face> &origMesh, std::vector<point3d> &destVert, std::vector<face> &destMesh, std::vector<vec3d> &destNorm) const;


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
    idxtype getNewVertex(const edge &e, std::vector<point3d> &vertList, const std::vector<face> &mesh, EdgeList &newVertList) const;



    [[deprecated]]
    void applyLoop(const face &t, const std::vector<point3d> &orig, std::vector<size_t> &valence, std::vector<point3d> &dest) const;

};