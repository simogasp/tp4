/**
 * It scales the model to unitary size by translating it to the origin and
 * scaling it to fit in a unit cube around the origin.
 * 
 * @return the scale factor used to transform the model
 */
float ObjModel::unitizeModel( )
{
	if ( !_vertices.empty( ) && !_mesh.empty( ) )
	{
		//****************************************
		// calculate model width, height, and 
		// depth using the bounding box
		//****************************************
		float w, h, d;
		w = fabs( _bb.pmax.x - _bb.pmin.x );
		h = fabs( _bb.pmax.y - _bb.pmin.y );
		d = fabs( _bb.pmax.z - _bb.pmin.z );

		cout << "size: w: " << w << " h " << h << " d " << d << endl;
		//****************************************
		// calculate center of the bounding box of the model 
		//****************************************
		point3d c = (_bb.pmax + _bb.pmin) * 0.5;

		//****************************************
		// calculate the unitizing scale factor as the 
		// maximum of the 3 dimensions
		//****************************************
		float scale = 2.0 / std::max( std::max( w, h ), d );

		cout << "scale: " << scale << " cx " << c.x << " cy " << c.y << " cz " << c.z << endl;

		// translate each vertex wrt to the center and then apply the scaling to the coordinate
		for ( int i = 0; i < _vertices.size( ); i++ )
		{
			//****************************************
			// translate the vertex
			//****************************************
			_vertices[i].translate( -c.x, -c.y, -c.z );

			//****************************************
			// apply the scaling
			//****************************************
			_vertices[i].scale( scale );

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

void ObjModel::release( ) { }

/**
 * It parses a line of the OBJ file containing a face and it return the result. 
 * NB: it only recover the indices, it discard normal and texture indices
 * 
 * @param toParse the string to parse in the OBJ format for a face (f v/vt/vn v/vt/vn v/vt/vn) and its variants
 * @param out the 3 indices for the face
 * @return true if the parse was successful
 */
bool ObjModel::parseFaceString( const string &toParse, face &out ) const
{
	if ( toParse.c_str( )[0] == 'f' )
	{
		idxtype a;
		// now check the different formats: %d, %d//%d, %d/%d, %d/%d/%d
		if ( strstr( toParse.c_str( ), "//" ) )
		{
			// v//n
			return ( sscanf( toParse.c_str( ), "f %u//%u %u//%u %u//%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a ) == 6);
		}
		else if ( sscanf( toParse.c_str( ), "f %u/%u/%u", &a, &a, &a ) == 3 )
		{
			// v/t/n
			return ( sscanf( toParse.c_str( ), "f %u/%u/%u %u/%u/%u %u/%u/%u", &(out.v1), &a, &a, &(out.v2), &a, &a, &(out.v3), &a, &a ) == 9);
		}
		else if ( sscanf( toParse.c_str( ), "f %u/%u", &a, &a ) == 2 )
		{
			// v/t .
			return ( sscanf( toParse.c_str( ), "f %u/%u %u/%u %u/%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a ) == 6);
		}
		else
		{
			// v
			sscanf( toParse.c_str( ), "f %u %u %u", &(out.v1), &(out.v2), &(out.v3) );
			PRINTVAR( out );
			return ( sscanf( toParse.c_str( ), "f %u %u %u", &(out.v1), &(out.v2), &(out.v3) ) == 3);
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

void ObjModel::flatDraw( ) const
{
	glShadeModel( GL_SMOOTH );

	// for each triangle draw the vertices and the normals
	for ( int i = 0; i < _mesh.size( ); i++ )
	{
		glBegin( GL_TRIANGLES );
		//compute the normal of the triangle
		vec3d n;
		computeNormal( _vertices[_mesh[i].v1], _vertices[(int) _mesh[i].v2], _vertices[(int) _mesh[i].v3], n );
		glNormal3fv( (float*) &n );

		glVertex3fv( (float*) &_vertices[_mesh[i].v1] );

		glVertex3fv( (float*) &_vertices[_mesh[i].v2] );

		glVertex3fv( (float*) &_vertices[_mesh[i].v3] );

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
	glDrawElements( GL_TRIANGLES, _mesh.size( ) * VERTICES_PER_TRIANGLE, GL_UNSIGNED_INT, (idxtype*) & _mesh[0] );

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

	glDrawElements( GL_TRIANGLES, _subMesh.size( ) * VERTICES_PER_TRIANGLE, GL_UNSIGNED_SHORT, (idxtype*) & _subMesh[0] );


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
	//	PRINTVAR(valence[t.v1]);

	valence[t.v2]++;
	destVert[t.v2] += (0.625f * origVert[t.v2] + 0.1875f * origVert[t.v1] + 0.1875f * origVert[t.v3]);

	valence[t.v3]++;
	destVert[t.v3] += (0.625f * origVert[t.v3] + 0.1875f * origVert[t.v2] + 0.1875f * origVert[t.v1]);
}
