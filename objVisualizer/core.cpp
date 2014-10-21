/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "core.hpp"

#include <cmath>
#include <algorithm>

void v3f::normalize()
{
    float n = norm();
	if (n > std::numeric_limits<float>::epsilon() )
	{
		x /= n;
		y /= n;
		z /= n;
	}
}

float v3f::dot(const v3f &v) const
{
    return (x*v.x + y*v.y + z*v.z);
}

float v3f::norm() const
{
    return sqrt(x*x+ y*y + z*z);
}

/**
 * Translate the vector
 * @param x the delta x of the translation
 * @param y the delta y of the translation
 * @param z the delta z of the translation
 */
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

void v3f::min( const v3f& a )
{
	x = std::min( x, a.x );
	y = std::min( y, a.y );
	z = std::min( z, a.z );
}

void v3f::max( const v3f& a )
{
	x = std::max( x, a.x );
	y = std::max( y, a.y );
	z = std::max( z, a.z );
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



//**************************** tindex ***********************//

// ********** SUM
tindex tindex::operator +(const tindex& a) const
{
	return tindex(v1 + a.v1, v2 + a.v2, v3 + a.v3);

}
tindex& tindex::operator +=(const tindex& a)
{
	v1 += a.v1;
	v2 += a.v2;
	v3 += a.v3;
	return *this;
}

tindex tindex::operator +(const GLushort a) const
{
	return tindex(v1 + a, v2 + a, v3 + a);
}
tindex& tindex::operator +=(const GLushort a)
{
	v1 += a;
	v2 += a;
	v3 += a;
	return *this;
}

// ********** DIFFERENCE
tindex tindex::operator -(const tindex& a) const
{
	return tindex(v1 - a.v1, v2 - a.v2, v3 - a.v3);
}
tindex& tindex::operator -=(const tindex& a)
{
	v1 -= a.v1;
	v2 -= a.v2;
	v3 -= a.v3;
	return *this;
}

tindex tindex::operator -(const GLushort a) const
{
	return tindex(v1 - a, v2 - a, v3 - a);
}
tindex& tindex::operator -=(const GLushort a)
{
	v1 -= a;
	v2 -= a;
	v3 -= a;
	return *this;
}

// ********** MULTIPLICATION
tindex tindex::operator *(const tindex& a) const
{
	return tindex(v1 * a.v1, v2 * a.v2, v3 * a.v3);
}
tindex& tindex::operator *=(const tindex& a)
{
	v1 *= a.v1;
	v2 *= a.v2;
	v3 *= a.v3;
	return *this;
}

tindex tindex::operator *(const GLushort a) const
{
	return tindex(v1 * a, v2 * a, v3 * a);
}
tindex& tindex::operator *=(const GLushort a)
{
	v1 *= a;
	v2 *= a;
	v3 *= a;
	return *this;
}

bool tindex::operator==( const tindex& r) const
{
	return ( (v1 == r.v1) && (v2 == r.v2) && (v3 == r.v3));
}
