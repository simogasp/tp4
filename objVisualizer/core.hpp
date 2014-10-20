/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _CORE_HPP_
#define	_CORE_HPP_

// for mac osx
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
// only for windows
#ifdef _WIN32
#include <windows.h>
#endif
// for windows and linux
#include <GL/gl.h>
#endif

#include <vector>
#include <map>
#include <iostream>
#include <functional>

#define DEBUGGING 1

#if DEBUGGING
#define PRINTVAR( a ) std::cout << #a << " = " << a << std::endl << std::endl;
#else
#define PRINTVAR( a )
#endif

typedef std::pair<GLushort,GLushort> edge;

inline GLushort sum( edge e)
{
    return (e.first + e.second);
}

inline GLushort min( edge e)
{
    return ((e.first > e.second) ? (e.second) : (e.first));
}


struct edgeCompare
{
  bool operator() (const edge &a, const edge &b) const
  {
//      return !( ( (a.first == b.first) && (a.second == b.second) ) ||
//               ( (a.first == b.second) && (a.second == b.first) ) );

      // two edges are equal if their sum is equal and their min is the same
      // our ordering policy is that first is the sum to be compared and then if the sum is
      // equal we use the ordering on the min. This guarantees that when edgeCompare is used
      // in map to detect equal keys it really detects the same edge independently of the index
      // order. edgeCompare is indeed used with reflexively: 2 edges are equivalent if !comp(a,b) && !comp(b,a)
      return ( (sum(a) < sum(b)) ||
               ((sum(a) == sum(b)) && (min(a) < min(b))));
  }
};

// to be used with unordered
//struct edgeHash
//{
//  size_t operator() (const edge &a ) const
//  {
//      std::hash<std::string> fun;
//      return (fun( (a.first > a.second) ? (std::to_string(a.second)+std::to_string(a.first) ) :
//                                          (std::to_string(a.first)+std::to_string(a.second)) ));
//  }
//};

typedef std::map< edge, GLushort, edgeCompare > _edgeList;

inline std::ostream& operator<<(std::ostream& os, const edge& p)
{
    return os << "["<< p.first << "," << p.second  <<"]";
}

inline std::ostream& operator<<(std::ostream& os, const _edgeList& l)
{
    os << std::endl;
    for( _edgeList::const_iterator it = l.begin(); it != l.end(); ++it )
        os << "\t" << it->first << "\t" << it->second << std::endl;
    return os;
}


class EdgeList
{

public:
    EdgeList() {}

    void add( edge e, GLushort idx )
    {
        list[e]=idx;
    }

    bool contains(edge e)
    {
        return (list.find( e ) != list.end());
    }

    GLushort getIndex(edge e)
    {
        return (list[e]);
    }

	friend std::ostream& operator<<(std::ostream& os, const EdgeList& l);


private:
    _edgeList list;

};

inline std::ostream& operator<<(std::ostream& os, const EdgeList& l)
{
	return (os << l.list);
}


struct v3f
{
    float x;
    float y;
    float z;

    v3f(float x, float y, float z) : x(x), y(y), z(z) {}

    v3f() : x(0), y(0), z(0) {}

    v3f(float a[3]) : x(a[0]), y(a[1]), z(a[2]) {}

    void normalize();

    float dot(const v3f &v);

    float norm() const;

    void translate( float x, float y, float z );

    void translate( v3f t );

    void scale( v3f t );

    void scale( float x, float y, float z );

    void scale( float a );

	void min( const v3f& a );

	void max( const v3f& a );

    v3f cross(const v3f& v);

    v3f cross(const float v[3]);

    v3f operator +(const v3f& a) const;
    v3f& operator +=(const v3f& a);

    v3f operator +(const float a[3]) const;
    v3f& operator +=(const float a[3]);

    v3f operator +(const float a) const;
    v3f& operator +=(const float a);

    // subtraction

    v3f operator -(const v3f& a) const;
    v3f& operator -=(const v3f& a);

    v3f operator -(const float a[3]) const;
    v3f& operator -=(const float a[3]);

    v3f operator -(const float a) const;
    v3f& operator -=(const float a);

    // element-wise product

    v3f operator *(const v3f& a) const;
    v3f& operator *=(const v3f& a);

    v3f operator *(const float a[3]) const;
    v3f& operator *=(const float a[3]);

    v3f operator *(const float a) const;
    v3f& operator *=(const float a);




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

	tindex() : v1(0), v2(0), v3(0) {}

    tindex(GLushort v1, GLushort v2, GLushort v3) : v1(v1), v2(v2), v3(v3) {}

	tindex operator +(const tindex& a) const;
	tindex& operator +=(const tindex& a);

	tindex operator +(const GLushort a) const;
	tindex& operator +=(const GLushort a);

	tindex operator -(const tindex& a) const;
	tindex& operator -=(const tindex& a);

	tindex operator -(const GLushort a) const;
	tindex& operator -=(const GLushort a);

	tindex operator *(const tindex& a) const;
	tindex& operator *=(const tindex& a);

	tindex operator *(const GLushort a) const;
	tindex& operator *=(const GLushort a);

	bool operator==(const tindex& rhs) const;
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

#endif //_CORE_HPP_
