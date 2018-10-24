#pragma once

#include "vec.h"
#include "aabb.h"
#include "particle.h"

struct node
{
    node* next[8];
    vec pos;
    real mass;
    particle *p;
    AABB aabb;
    real size;
    bool isLeaf;

    node()
    {
        reset();
    }

    void reset()
    {
        for(int i = 0; i < 8; ++i) {
            next[i] = 0;
        }
        pos.setZero();
        mass = real(0);
        p = 0;
        aabb.reset();
        size = real(0);
        isLeaf = true;
    }
};
