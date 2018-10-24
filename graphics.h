#ifndef GRAPHICS_INCLUDE
#define GRAPHICS_INCLUDE

#include "aabb.h"
#include "node.h"
#include "particle.h"

void initPaint(int w, int h, AABB& aabb);
void paint(const particle *particles, const int n, node * const root);
void quitPaint();

#endif
