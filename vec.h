#pragma once

#include "bhtmath.h"

struct vec
{
    union
    {
        real data[3];

        struct
        {
            real x;
            real y;
            real z;
        };
    };

    void set(real x_, real y_, real z_)
    {
        x = x_;
        y = y_;
        z = z_;
    }

    void setZero()
    {
        set(0, 0, 0);
    }

    vec()
    {
        setZero();
    }

    vec(const vec& other)
    {
        set(other.x, other.y, other.z);
    }

    vec(real x_, real y_, real z_)
    {
        set(x_, y_, z_);
    }

    vec& operator=(const vec &other)
    {
        set(other.x, other.y, other.z);
        return (*this);
    }

    vec operator+(const vec &other) const
    {
        return vec(x + other.x, y + other.y, z + other.z);
    }

    vec operator-(const vec &other) const 
    {
        return vec(x - other.x, y - other.y, z - other.z);
    }

    vec operator*(real scaling) const
    {
        return vec(x * scaling, y * scaling, z * scaling);
    }

    vec operator/(real scaling) const
    {
        const real inv = real(1) / scaling;
        return vec(x * inv, y * inv, z * inv);
    }

    vec operator+(real val) const
    {
        return vec(x + val, y + val, z + val);
    }

    vec operator-(real val) const
    {
        return vec(x - val, y - val, z - val);
    }

    real operator*(const vec &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    vec operator%(const vec &other) const
    {
        return vec(
            y*other.z - z*other.y,
            z*other.x - x*other.z,
            x*other.y - y*other.x
        );
    }

    real sqDist(const vec& other) const
    {
        const vec diff = (*this) - other;
        return diff * diff;
    }
    
    real dist(const vec& other) const
    {
        const real squareDistance = sqDist(other);
        return Math::sqrt(squareDistance);
    }

    real length() const
    {
        return Math::sqrt(x*x+y*y+z*z);
    }

    real normalize()
    {
        real len = length();
        (*this) = (*this) / len;
        return len;
    }

    real getNormalized(vec &out) const
    {
        out = (*this);
        return out.normalize();
    }
};
