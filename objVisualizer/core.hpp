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
#include <iostream>

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

    v3f cross(const v3f& v);

    v3f cross(const float v[3]);


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

#endif //_CORE_HPP_