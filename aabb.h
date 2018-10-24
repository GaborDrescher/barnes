#pragma once
#include "vec.h"
#include "bhtmath.h"

struct AABB
{
    vec min;
    vec max;

    void reset()
    {
        real minval = REAL_MIN;
        real maxval = REAL_MAX;

        min.set(maxval, maxval, maxval);
        max.set(minval, minval, minval);
    }

    AABB()
    {
        reset();
    }

    void addVec(const vec& v)
    {
        for(int i = 0; i < 3; ++i) {
            min.data[i] = Math::min(min.data[i], v.data[i]);
            max.data[i] = Math::max(max.data[i], v.data[i]);
        }
    }

    bool isInside(const vec& v) const
    {
        for(int i = 0; i < 3; ++i) {
            if(v.data[i] < min.data[i] || v.data[i] > max.data[i]) {
                return false;
            }
        }

        return true;
    }

    bool isInside(const vec& v, const real scale) const
    {
        for(int i = 0; i < 3; ++i) {
            if(v.data[i] < min.data[i] * scale || v.data[i] > max.data[i] * scale) {
                return false;
            }
        }

        return true;
    }

    real getLongestSide() const
    {
        real longest = REAL_MIN;
        for(int i = 0; i < 3; ++i) {
            longest = Math::max(longest, max.data[i] - min.data[i]);
        }

        return longest;
    }

    vec getCenter() const
    {
        return (min + max) * real(0.5);
    }
};
