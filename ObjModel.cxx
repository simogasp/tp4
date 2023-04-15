/**
 * @file ObjModel.cxx
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

#include <algorithm>
#include <regex>
#include <array>

// lambda to convert the string matches to a face
auto face_from_match(const std::string& vst1, const std::string& vst2, const std::string& vst3) -> face
{
    const auto v1 = static_cast<idxtype>(std::stoi(vst1));
    const auto v2 = static_cast<idxtype>(std::stoi(vst2));
    const auto v3 = static_cast<idxtype>(std::stoi(vst3));
    return {v1, v2, v3};
}

face parseFaceString(const std::string &toParse)
{
    const auto res = parseFaceStringRegex(toParse);
    if(!res.has_value())
        throw std::runtime_error("Error while reading line: " + toParse);

    return res.value();
}

std::optional<face> parseFaceStringRegex( const std::string &toParse)
{
    static const std::array<std::regex, 4> face_regexes {
          // regex to match a face from a string in the format f v1 v2 v3
          std::regex(R"(f\s+(\d+)\s+(\d+)\s+(\d+))"),
          // regex to match a face from a string in the format f v1/t1 v2/t2 v3/t3
          std::regex(R"(f\s+(\d+)(?:\/\d+){1}\s+(\d+)(?:\/\d+){1}\s+(\d+)(?:\/\d+){1})"),
          // regex to match a face from a string in the format f v1/t1/n1 v2/t2/n2 v3/t3/n3
          std::regex(R"(f\s+(\d+)(?:\/\d+){2}\s+(\d+)(?:\/\d+){2}\s+(\d+)(?:\/\d+){2})"),
          // regex to match a face from a string in the format f v1//n1 v2//n2 v3//n3
          std::regex (R"(f\s+(\d+)(?:\/\/\d+){1}\s+(\d+)(?:\/\/\d+){1}\s+(\d+)(?:\/\/\d+){1})")
    };

    // early exit if the string is empty or does not start with 'f'
    if (toParse.empty() || toParse[0] != 'f' )
    {
        return {};
    }

    std::smatch face_match;

    // try to match the string to each of the regexes
    for(const auto& face_regex : face_regexes)
    {
        if(std::regex_search(toParse, face_match, face_regex))
        {
            return face_from_match(face_match[1], face_match[2], face_match[3]);
        }
    }

    // if no match was found, return an empty optional
    return {};
}

point3d parseVertexString(const std::string &toParse)
{
    const auto res = parseVertexStringRegex(toParse);
    if(!res.has_value())
        throw std::runtime_error("Error while reading line: " + toParse);

    return res.value();
}

std::optional<point3d> parseVertexStringRegex(const std::string &toParse)
{
    static const std::regex vertex_regex(R"(v\s+([+-]?(?:[0-9]*[.])?[0-9]+)\s+([+-]?(?:[0-9]*[.])?[0-9]+)\s+([+-]?(?:[0-9]*[.])?[0-9]+))");

    // early exit if the string is empty or does not start with 'v'
    if (toParse.empty() || toParse[0] != 'v' )
    {
        return {};
    }

    std::smatch vertex_match;

    // try to match the string to the regex
    if(std::regex_search(toParse, vertex_match, vertex_regex))
    {
        const auto x = std::stof(vertex_match[1]);
        const auto y = std::stof(vertex_match[2]);
        const auto z = std::stof(vertex_match[3]);
        return point3d(x, y, z);
    }

    // if no match was found, return an empty optional
    return {};
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
* Render the model according to the provided parameters
* @param params The rendering parameters
*/
void ObjModel::render( const RenderingParameters &params )
{
    // if we need to draw the original model
    if ( !params.subdivision )
    {
        // draw it
        draw( _vertices, _mesh, _normals, params );
        // draw the normals
        if ( params.normals )
        {
            drawNormals( _vertices, _normals );
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
            vector<point3d> tmpVert;        //!< a temporary list of vertices used in the iterations
            vector<face> tmpMesh;           //!< a temporary mesh used in the iterations

            if(( _currentSubdivLevel == 0 ) || ( _currentSubdivLevel > params.subdivLevel ) )
            {
                // start from the beginning
                _currentSubdivLevel = 0;
                tmpVert = _vertices;
                tmpMesh = _mesh;
            }
            else
            {
                // start from the current level
                tmpVert = _subVert;
                tmpMesh = _subMesh;
            }

            // apply the proper subdivision iterations
            for( ; _currentSubdivLevel < params.subdivLevel; ++_currentSubdivLevel)
            {
                cerr << "[Loop subdivision] iteration " << _currentSubdivLevel << endl;
                loopSubdivision( tmpVert, tmpMesh, _subVert, _subMesh, _subNorm );
                // swap unless it's the last iteration
                if( _currentSubdivLevel < ( params.subdivLevel - 1) )
                {
                    tmpVert = _subVert;
                    tmpMesh = _subMesh;
                }
            }
        }

        draw( _subVert, _subMesh, _subNorm, params );
        if ( params.normals )
        {
            drawNormals( _subVert, _subNorm );
        }
    }
}

/**
 * Draw the model
 * 
 * @param vertices list of vertices
 * @param indices list of faces
 * @param vertexNormals list of normals
 * @param params Rendering parameters
 */
void ObjModel::draw( const std::vector<point3d> &vertices, const std::vector<face> &indices, std::vector<vec3d> &vertexNormals, const RenderingParameters &params ) const
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

void ObjModel::drawSolid( const vector<point3d> &vertices, const vector<face> &indices, vector<vec3d> &vertexNormals, const RenderingParameters &params ) const
{
    if ( params.useIndexRendering )
    {
        drawSmoothFaces( vertices, indices, vertexNormals, params );
    }
    else
    {
        drawFlatFaces( vertices, indices, params );
    }
}

/**
 * Draw the normals at each vertex
 * @param vertices The list of vertices 
 * @param vertexNormals The list of associated normals
 */
void ObjModel::drawNormals( const std::vector<point3d> &vertices, const std::vector<vec3d>& vertexNormals ) const
{
    glDisable( GL_LIGHTING );

    glColor3f( 0.8f, .0f, .0f );
    glLineWidth( 2 );

    for(std::size_t i = 0; i < vertices.size(); ++i)
    {
        glBegin( GL_LINES );

        const auto v = vertices[i];
        const auto n = vertexNormals[i];

        vec3d newP = v +  0.05f * n;
        glVertex3fv( (float*) &v );

        glVertex3f( newP.x, newP.y, newP.z );

        glEnd( );
    }
    glEnable( GL_LIGHTING );
}

/**
 * It scales the model to unitary size by translating it to the origin and
 * scaling it to fit in a unit cube around the origin.
 * 
 * @return the scale factor used to transform the model
 */
float ObjModel::unitizeModel( )
{
    if ( _vertices.empty( ) || _mesh.empty( ) )
    {
        return .0f;
    }

    //****************************************
    // calculate model width, height, and
    // depth using the bounding box
    //****************************************
    const float w = fabs( _bb.pmax.x - _bb.pmin.x );
    const float h = fabs( _bb.pmax.y - _bb.pmin.y );
    const float d = fabs( _bb.pmax.z - _bb.pmin.z );

    cout << "size: w: " << w << " h " << h << " d " << d << endl;
    //****************************************
    // calculate center of the bounding box of the model
    //****************************************
    const point3d c = (_bb.pmax + _bb.pmin) * 0.5;

    //****************************************
    // calculate the unitizing scale factor as the
    // maximum of the 3 dimensions
    //****************************************
    const auto scale = static_cast<float>(2.f / std::max(std::max(w, h), d));

    cout << "scale: " << scale << " cx " << c.x << " cy " << c.y << " cz " << c.z << endl;

    // translate each vertex wrt to the center and then apply the scaling to the coordinate
    for(auto& v : _vertices)
    {
        //****************************************
        // translate the vertex
        //****************************************
        v.translate( -c.x, -c.y, -c.z );

        //****************************************
        // apply the scaling
        //****************************************
        v.scale( scale );

    }


    //****************************************
    // update the bounding box, ie translate and scale the 6 coordinates
    //****************************************
    _bb.pmax = (_bb.pmax - c) * scale;
    _bb.pmin = (_bb.pmin - c) * scale;


    cout << "New bounding box : pmax=" << _bb.pmax << "  pmin=" << _bb.pmin << endl;

    return scale;

}


//*****************************************************************************
//*                        DEPRECATED FUNCTIONS
//*****************************************************************************

// to be deprecated

void ObjModel::flatDraw( ) const
{
    glShadeModel( GL_SMOOTH );

    // for each triangle draw the vertices and the normals
    for(const auto &face : _mesh)
    {
        glBegin( GL_TRIANGLES );
        //compute the normal of the triangle
        vec3d n;
        computeNormal( _vertices[face.v1], _vertices[face.v2], _vertices[face.v3], n );
        glNormal3fv( (float*) &n );

        glVertex3fv( (float*) &_vertices[face.v1] );

        glVertex3fv( (float*) &_vertices[face.v2] );

        glVertex3fv( (float*) &_vertices[face.v3] );

        glEnd( );
    }

}

// to be deprecated

void ObjModel::drawWireframe( ) const
{

    drawWireframe( _vertices, _mesh, RenderingParameters( ) );

}

// to be deprecated

void ObjModel::indexDraw( ) const
{
    glShadeModel( GL_SMOOTH );
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
    glNormalPointer( GL_FLOAT, 0, (float*) &_normals[0] );

    //****************************************
    // Index pointer to normal array
    //****************************************
    glVertexPointer( COORD_PER_VERTEX, GL_FLOAT, 0, (float*) &_vertices[0] );

    //****************************************
    // Draw the triangles
    //****************************************
    glDrawElements( GL_TRIANGLES, static_cast<GLsizei>(_mesh.size( )) * VERTICES_PER_TRIANGLE, GL_UNSIGNED_INT, (idxtype*) & _mesh[0] );

    //****************************************
    // Disable vertex arrays
    //****************************************
    glDisableClientState( GL_VERTEX_ARRAY ); // disable vertex arrays

    //****************************************
    // Disable normal arrays
    //****************************************
    glDisableClientState( GL_NORMAL_ARRAY );
}

// to be deprecated

void ObjModel::drawSubdivision( )
{
    if ( _subMesh.empty( ) || _subNorm.empty( ) || _subVert.empty( ) )
    {
        loopSubdivision( _vertices, _mesh, _subVert, _subMesh, _subNorm );
    }

    glShadeModel( GL_SMOOTH );

    glEnableClientState( GL_NORMAL_ARRAY );
    glEnableClientState( GL_VERTEX_ARRAY );

    glNormalPointer( GL_FLOAT, 0, (float*) &_subNorm[0] );
    glVertexPointer( COORD_PER_VERTEX, GL_FLOAT, 0, (float*) &_subVert[0] );

    glDrawElements( GL_TRIANGLES, static_cast<GLsizei>(_subMesh.size( )) * VERTICES_PER_TRIANGLE, GL_UNSIGNED_SHORT, (idxtype*) & _subMesh[0] );


    glDisableClientState( GL_VERTEX_ARRAY ); // disable vertex arrays
    glDisableClientState( GL_NORMAL_ARRAY );

    drawWireframe( _subVert, _subMesh, RenderingParameters( ) );

}

// to be removed
void ObjModel::applyLoop( const face &t, const std::vector<point3d> &origVert, std::vector<size_t> &valence, std::vector<point3d> &destVert ) const
{
    // 5/8 V + 3/8 sum(V_i))
    // in this case since we are summing each face the other vertices are counted
    // twice, so we use 3/16 instead of 3/8
    valence[t.v1]++;
    destVert[t.v1] += (0.625f * origVert[t.v1] + 0.1875f * origVert[t.v2] + 0.1875f * origVert[t.v3]);
    //    PRINTVAR(valence[t.v1]);

    valence[t.v2]++;
    destVert[t.v2] += (0.625f * origVert[t.v2] + 0.1875f * origVert[t.v1] + 0.1875f * origVert[t.v3]);

    valence[t.v3]++;
    destVert[t.v3] += (0.625f * origVert[t.v3] + 0.1875f * origVert[t.v2] + 0.1875f * origVert[t.v1]);
}
