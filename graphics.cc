#include "graphics.h"
#include "vec.h"
#include <SDL.h>
#include <SDL_opengl.h>

static real boxlen;

void initPaint(int w, int h, AABB& aabb)
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_SetVideoMode(w, h, 0, SDL_OPENGL | SDL_HWSURFACE);

    //real ratio = (real) width / (real) height;
    glShadeModel(GL_SMOOTH);
    glClearColor(0, 0, 0, 0);
    glEnable(GL_DEPTH_TEST);
    
    //glEnable(GL_POINT_SMOOTH);
    //glEnable(GL_LINE_SMOOTH);
    glPointSize(1);

    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, w/(real)h, 1, 100);//near/far plane

    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
    gluLookAt(3, 3, -5, 0, 0, 0, 0, 1, 0);

    boxlen = real(0);
    boxlen = Math::max(boxlen, aabb.max.data[0] - aabb.min.data[0]);
    boxlen = Math::max(boxlen, aabb.max.data[1] - aabb.min.data[1]);
    boxlen = Math::max(boxlen, aabb.max.data[2] - aabb.min.data[2]);
}

static void drawTreeRek(node * const root, int depth)
{


    if(root == 0) return;
    if(root->isLeaf && root->p == 0) return;

    if(depth < 5) {
        AABB &aabb = root->aabb;
        real invBoxlen = real(1) / boxlen;

        //vec center = (aabb.min + aabb.max) * real(0.5);

        glColor3f(0,1,0);
        glBegin(GL_QUAD_STRIP);
        glVertex3f(aabb.max.data[0] * invBoxlen, aabb.max.data[1] * invBoxlen, aabb.max.data[2] * invBoxlen);
        glVertex3f(aabb.max.data[0] * invBoxlen, aabb.max.data[1] * invBoxlen, aabb.min.data[2] * invBoxlen);
        glVertex3f(aabb.min.data[0] * invBoxlen, aabb.max.data[1] * invBoxlen, aabb.max.data[2] * invBoxlen);
        glVertex3f(aabb.min.data[0] * invBoxlen, aabb.max.data[1] * invBoxlen, aabb.min.data[2] * invBoxlen);
        glVertex3f(aabb.min.data[0] * invBoxlen, aabb.min.data[1] * invBoxlen, aabb.max.data[2] * invBoxlen);
        glVertex3f(aabb.min.data[0] * invBoxlen, aabb.min.data[1] * invBoxlen, aabb.min.data[2] * invBoxlen);
        glVertex3f(aabb.max.data[0] * invBoxlen, aabb.min.data[1] * invBoxlen, aabb.max.data[2] * invBoxlen);
        glVertex3f(aabb.max.data[0] * invBoxlen, aabb.min.data[1] * invBoxlen, aabb.min.data[2] * invBoxlen);
        glVertex3f(aabb.max.data[0] * invBoxlen, aabb.max.data[1] * invBoxlen, aabb.max.data[2] * invBoxlen);
        glVertex3f(aabb.max.data[0] * invBoxlen, aabb.max.data[1] * invBoxlen, aabb.min.data[2] * invBoxlen);
        glEnd();

    }

    depth += 1;
    for(int i = 0; i < 8; ++i)
    {
        drawTreeRek(root->next[i], depth);
    }
}

static void drawTree(node * const head)
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawTreeRek(head, 0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

static bool doDrawTree = false;

void paint(const particle *particles, const int n, node * const root)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //paint here
    glColor3f(1,0,0);
    glBegin(GL_POINTS);
    real invBoxlen = real(1) / boxlen;
    for(int i = 0; i < n; ++i)
    {
        const particle& p = particles[i];
        real x = p.pos.data[0] * invBoxlen;
        real y = p.pos.data[1] * invBoxlen;
        real z = (p.pos.data[2] * invBoxlen);

        //std::cout << x << " " << y << " " << z << std::endl;
        
        glVertex3f(x,y,z);
    }
    glEnd();

    if(doDrawTree && root != 0) {
        drawTree(root);
    }
    
    SDL_GL_SwapBuffers();

    SDL_Event event;
    while(SDL_PollEvent(&event))
    {
        switch (event.type)
        {
            case SDL_QUIT:
                exit(0);
                break;
            case SDL_KEYDOWN:
                doDrawTree = !doDrawTree;
                break;
        }
    }
}

void quitPaint()
{
    SDL_Quit();
}

