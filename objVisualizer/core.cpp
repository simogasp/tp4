/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "core.hpp"

#include <cmath>

void v3f::normalize()
{
    float n = norm();
    x /= n;
    y /= n;
    z /= n;
}

float v3f::dot(const v3f &v)
{
    return (x*v.x + y*v.y + z*v.z);
}

float v3f::norm() const
{
    return sqrt(x*x+ y*y + z*z);
}

void v3f::translate( float x, float y, float z )
{
    this->x += x;
    this->y += y;
    this->z += z;
}

void v3f::translate( v3f t )
{
    x += t.x;
    y += t.y;
    z += t.z;
}

void v3f::scale( v3f t )
{
    x *= t.x;
    y *= t.y;
    z *= t.z;
}

void v3f::scale( float x, float y, float z )
{
    scale(v3f(x,y,z));
}

void v3f::scale( float a )
{
    scale(v3f(a,a,a));
}

v3f v3f::operator +(const v3f& a) const
{
    return v3f(x + a.x, y + a.y, z + a.z);
}

v3f& v3f::operator +=(const v3f& a)
{
    x += a.x;
    y += a.y;
    z += a.z;
    return *this;
}

v3f v3f::operator +(const float a[3]) const
{
    return v3f(x + a[0], y + a[1], z + a[2]);
}

v3f& v3f::operator +=(const float a[3])
{
    x += a[0];
    y += a[1];
    z += a[2];
    return *this;
}

v3f v3f::operator +(const float a) const
{
    return v3f(x + a, y + a, z + a);
}

v3f& v3f::operator +=(const float a)
{
    x += a;
    y += a;
    z += a;
    return *this;
}

// subtraction

v3f v3f::operator -(const v3f& a) const
{
    return v3f(x - a.x, y - a.y, z - a.z);
}

v3f& v3f::operator -=(const v3f& a)
{
    x -= a.x;
    y -= a.y;
    z -= a.z;
    return *this;
}

v3f v3f::operator -(const float a[3]) const
{
    return v3f(x - a[0], y - a[1], z - a[2]);
}

v3f& v3f::operator -=(const float a[3])
{
    x -= a[0];
    y -= a[1];
    z -= a[2];
    return *this;
}

v3f v3f::operator -(const float a) const
{
    return v3f(x - a, y - a, z - a);
}

v3f& v3f::operator -=(const float a)
{
    x -= a;
    y -= a;
    z -= a;
    return *this;
}

// element-wise product

v3f v3f::operator *(const v3f& a) const
{
    return v3f(x * a.x, y * a.y, z * a.z);
}

v3f& v3f::operator *=(const v3f& a)
{
    x *= a.x;
    y *= a.y;
    z *= a.z;
    return *this;
}

v3f v3f::operator *(const float a[3]) const
{
    return v3f(x * a[0], y * a[1], z * a[2]);
}

v3f& v3f::operator *=(const float a[3])
{
    x *= a[0];
    y *= a[1];
    z *= a[2];
    return *this;
}

v3f v3f::operator *(const float a) const
{
    return v3f(x * a, y * a, z * a);
}

v3f& v3f::operator *=(const float a)
{
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

v3f v3f::cross(const v3f& v)
{
    return v3f(y*v.z - z*v.y,
               z*v.x - x*v.z,
               x*v.y - y*v.x);
}

v3f v3f::cross(const float v[3])
{
    return v3f(y*v[2] - z*v[1],
               z*v[0] - x*v[2],
               x*v[1] - y*v[0]);
}
