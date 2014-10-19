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

// for mac osx
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
// only for windows
#ifdef _WIN32
#include <windows.h>
#endif
// for windows and linux
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/freeglut.h>
#endif


#include <vector>
#include <ostream>
#include <math.h>

#define VERTICES_PER_TRIANGLE 3
#define COORD_PER_VERTEX 3
#define TOTAL_FLOATS_IN_TRIANGLE (VERTICES_PER_TRIANGLE*COORD_PER_VERTEX)

#define DEBUGGING 1

#if DEBUGGING
#define PRINTVAR( a ) std::cout << #a << " = " << a << endl << endl;
#else
#define PRINTVAR( a )
#endif

 /**
  */
typedef struct BoundingBox
{
	float Xmax;
	float Ymax;
	float Zmax;
	float Xmin;
	float Ymin;
	float Zmin;
	
	BoundingBox( ) :
	Xmax( 0.0 ),
	Ymax( 0.0 ),
	Zmax( 0.0 ),
	Xmin( 0.0 ),
	Ymin( 0.0),
    Zmin( 0.0 ) {}
} BoundingBox;
 

struct v3f
{
    float x;
    float y;
    float z;

    v3f(float x, float y, float z) : x(x), y(y), z(z) {}

    v3f() : x(0), y(0), z(0) {}

    v3f(float a[3]) : x(a[0]), y(a[1]), z(a[2]) {}

    void normalize()
    {
        float n = norm();
        x /= n;
        y /= n;
        z /= n;
    }

    float dot(const v3f &v)
    {
        return (x*v.x + y*v.y + z*v.z);
    }

    float norm() const
    {
        return sqrt(x*x+ y*y + z*z);
    }

    void translate( float x, float y, float z )
    {
        this->x += x;
        this->y += y;
        this->z += z;
    }

    void translate( v3f t )
    {
        x += t.x;
        y += t.y;
        z += t.z;
    }

    void scale( v3f t )
    {
        x *= t.x;
        y *= t.y;
        z *= t.z;
    }

    void scale( float x, float y, float z )
    {
        scale(v3f(x,y,z));
    }

    void scale( float a )
    {
        scale(v3f(a,a,a));
    }

    v3f operator +(const v3f& a) const
    {
        return v3f(x + a.x, y + a.y, z + a.z);
    }

    v3f& operator +=(const v3f& a)
    {
        x += a.x;
        y += a.y;
        z += a.z;
        return *this;
    }

    v3f operator +(const float a[3]) const
    {
        return v3f(x + a[0], y + a[1], z + a[2]);
    }

    v3f& operator +=(const float a[3])
    {
        x += a[0];
        y += a[1];
        z += a[2];
        return *this;
    }

    v3f operator +(const float a) const
    {
        return v3f(x + a, y + a, z + a);
    }

    v3f& operator +=(const float a)
    {
        x += a;
        y += a;
        z += a;
        return *this;
    }

    // subtraction

    v3f operator -(const v3f& a) const
    {
        return v3f(x - a.x, y - a.y, z - a.z);
    }

    v3f& operator -=(const v3f& a)
    {
        x -= a.x;
        y -= a.y;
        z -= a.z;
        return *this;
    }

    v3f operator -(const float a[3]) const
    {
        return v3f(x - a[0], y - a[1], z - a[2]);
    }

    v3f& operator -=(const float a[3])
    {
        x -= a[0];
        y -= a[1];
        z -= a[2];
        return *this;
    }

    v3f operator -(const float a) const
    {
        return v3f(x - a, y - a, z - a);
    }

    v3f& operator -=(const float a)
    {
        x -= a;
        y -= a;
        z -= a;
        return *this;
    }

    // element-wise product

    v3f operator *(const v3f& a) const
    {
        return v3f(x * a.x, y * a.y, z * a.z);
    }

    v3f& operator *=(const v3f& a)
    {
        x *= a.x;
        y *= a.y;
        z *= a.z;
        return *this;
    }

    v3f operator *(const float a[3]) const
    {
        return v3f(x * a[0], y * a[1], z * a[2]);
    }

    v3f& operator *=(const float a[3])
    {
        x *= a[0];
        y *= a[1];
        z *= a[2];
        return *this;
    }

    v3f operator *(const float a) const
    {
        return v3f(x * a, y * a, z * a);
    }

    v3f& operator *=(const float a)
    {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }

    v3f cross(const v3f& v)
    {
        return v3f(y*v.z - z*v.y,
                   z*v.x - x*v.z,
                   x*v.y - y*v.x);
    }

    v3f cross(const float v[3])
    {
        return v3f(y*v[2] - z*v[1],
                   z*v[0] - x*v[2],
                   x*v[1] - y*v[0]);
    }


};

inline std::ostream& operator<<(std::ostream& os, const v3f& p)
{
    return os << "["<< p.x << "," << p.y << "," << p.z <<"]";
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<v3f>& p)
{
    os << std::endl;
    for(int i = 0; i<p.size(); ++i)
        os << "\t" << p[i] << std::endl;
    return os;
}

struct tindex
{
    GLushort v1;
    GLushort v2;
    GLushort v3;

    tindex(GLushort v1, GLushort v2, GLushort v3) : v1(v1), v2(v2), v3(v3) {}

};

inline std::ostream& operator<<(std::ostream& os, const tindex& p)
{
    return os << "["<< p.v1 << "," << p.v2 << "," << p.v3 <<"]";
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<tindex>& p)
{
    os << std::endl;
    for(int i = 0; i<p.size(); ++i)
        os << "\t" << p[i] << std::endl;
    return os;
}

typedef struct v3f point3d;
typedef struct v3f vec3d;
typedef struct tindex triangleIndex;

 
class ObjModel
{
  public: 
	ObjModel();			
		
		/**
		 * Calculate the normal of a triangular face defined by three points
		 * 
         * @param[in] coord1 the first vertex
         * @param[in] coord2 the second vertex
         * @param[in] coord3 the third vertex
         * @param[out] norm the normal
         */
        void computeNormal(const float coord1[3], const float coord2[3], const float coord3[3], float norm[3] );

        void computeNormal( const point3d& v1, const point3d& v2, const point3d& v3, vec3d &norm  ) const;

        float angleAtVertex( const point3d& v1, const point3d& v2, const point3d& v3 ) const;
		
		/**
		 *  Loads the model from file
		 * @param filename the OBJ file
		 * @return 
		 */
		int load(char *filename);	
		
		/**
		 * Draws the model with the opengl primitives
		 */
        void draw();

        void indexDraw();

        void flatDraw();

        void wireframetDraw();
	
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
		float* _normals;			// Stores the normals
		float* _triangles;			// Stores the triangles
		float* _vertices;			// Stores the points which make the object
        std::vector<triangleIndex> _indices;	// Stores the vertex indices for the triangles
        std::vector<point3d> _v;	// Stores the vertices
        std::vector<vec3d> _nt;      // Stores the normals for the triangles
        std::vector<vec3d> _nv;      // Stores the normals for the triangles
					
		long _numVertices;          // the actual number of loaded vertices
		long _numTriangles;         // the actual number of loaded faces
		

  private:
  
  
	/**
	 * Perform a first scan of the file in order to get the number of vertices and the number of faces
	 * @param filename
	 * @param vertexNum
	 * @param faceNum
	 */
	void firstScan(char* filename, long &vertexNum, long &faceNum);
	
  private:
	  
	// contains the bounding box of the model
	BoundingBox _bb;
 
};
 
#endif

