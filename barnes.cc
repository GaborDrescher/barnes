#include "config.h"
#include "aabb.h"
#include "bhtmath.h"
#include "graphics.h"
#include "node.h"
#include "nodealloc.h"
#include "particle.h"
#include "vec.h"
#include "xsrandom.h"

static void addForceSingle(vec& force, const vec& pos1, const vec& pos2, real m1, real m2)
{
    (void)m1;
    const vec r = pos2 - pos1;
    const real dist = Math::sqrt((r * r) + NEAR_EPS);
    const real invdist = real(1) / dist;
    const real invDistCube = invdist * invdist * invdist;

    force = force + (r * (m2 * invDistCube));
}

#ifdef USE_BRUTE_FORCE_ALGO
static void addForce(particle& a, particle& b)
{
    addForceSingle(a.force, a.pos, b.pos, a.mass, b.mass);
    addForceSingle(b.force, b.pos, a.pos, b.mass, a.mass);
}
#endif

static void integrate(particle& p, const real dt)
{
    vec tmp = p.force * (dt * GRAV_CONST);
    p.vel = p.vel + tmp;
    p.pos = p.pos + (p.vel * dt);

    p.force.setZero();
}

static void makeCube(AABB &aabb)
{
    vec center = aabb.getCenter();
    real longest = aabb.getLongestSide() * real(1.1);
    real half = longest * real(0.5);
    aabb.min = center - half;
    aabb.max = center + half;
}

static void createSubnodes(node *n)
{
    real size = n->size * real(0.5);
    AABB &aabb = n->aabb;
    vec mid = aabb.getCenter();
    vec vec1[2];
    vec vec2[2];
    vec1[0] = aabb.max;
    vec1[1] = mid;
    vec2[0] = mid;
    vec2[1] = aabb.min;

    for(int i = 0; i < 8; ++i) {
        node *current = getNewNode();
        vec &cmax = current->aabb.max;
        vec &cmin = current->aabb.min;

        for(int k = 0; k < 3; k++) {
            cmax.data[k] = vec1[((i >> k) & 1)].data[k];
            cmin.data[k] = vec2[((i >> k) & 1)].data[k];
        }

        current->size = size;
        n->next[i] = current;
    }
}

static node* chooseRight(node* root, particle *p)
{
    for(int i = 0; i < 8; ++i) {
        if(root->next[i]->aabb.isInside(p->pos)) {
            return root->next[i];
        }
    }

    return 0;//should not happen
}

static void addRek(node* root, particle *p)
{
    if(!root->isLeaf) {
        addRek(chooseRight(root, p), p);
    }
    else if(root->isLeaf && root->p != 0) {
        particle *old = root->p;
        root->isLeaf = false;
        root->p = 0;

        createSubnodes(root);
        addRek(chooseRight(root, old), old);
        addRek(chooseRight(root, p), p);
    }
    else /*if(root.isLeaf && root.p == null)*/{
        root->p = p;
    }
}

static void adjustTree(node *root)
{
    if(root->isLeaf && root->p == 0) {
        return;
    }
    
    if(root->p != 0)
    {
        root->mass = root->p->mass;
        root->pos = root->p->pos;
        return;
    }
    
    for(int i = 0; i < 8; ++i) {
        adjustTree(root->next[i]);
    }
    
    for(int i = 0; i < 8; ++i) {
        root->mass += root->next[i]->mass;
        vec tmp = root->next[i]->pos * root->next[i]->mass;
        root->pos = root->pos + tmp;
    }
    
    root->pos = root->pos / root->mass;
}

static node* addToTree(particle *particles, int n)
{
    node *root = getNewNode();

    for(int i = 0; i < n; ++i) {
        particle &p = particles[i];
        root->aabb.addVec(p.pos);
    }

    makeCube(root->aabb);
    root->size = root->aabb.max.data[0] - root->aabb.min.data[0];

    for(int i = 0; i < n; ++i) {
        addRek(root, particles + i);
    }

    adjustTree(root);
    
    return root;
}

static void calcTreeForce(node *root, particle& p, real theta)
{
    if(root->p != 0) {
        if(root->p == (&p)) return;
        addForceSingle(p.force, p.pos, root->pos, p.mass, root->mass);
    }
    else if(!root->isLeaf) {
        real s = root->size;
        real d = p.pos.dist(root->pos);
        real quot = s / d;
        if(quot < theta) {
            addForceSingle(p.force, p.pos, root->pos, p.mass, root->mass);
        }
        else {
            for(int i = 0; i < 8; ++i) {
                calcTreeForce(root->next[i], p, theta);
            }
        }
    }
}

static void destroyTree(node *head)
{
    if(head == 0) return;
    if(head->isLeaf) {
        destroyNode(head);
    }
    else {
        for(int i = 0; i < 8; ++i) {
            destroyTree(head->next[i]);
        }
    }
}

static real rndReal(XSRandom &rnd)
{
    uint64_t full = rnd.next();
    uint32_t ival = (uint32_t)(full ^ (full >> 32));
    double dval = ival / (double)UINT32_MAX;
    return (real)dval;
}

static void createScene(particle *particles, const int n, XSRandom &rnd)
{
    real centermass(2.0E5);
    real diam(9.461E2);

    for(int i = 0; i < n; ++i)
    {
        particle &p = particles[i];
        p.mass = real(2.0E11);

        real x = rndReal(rnd) * real(2.0) - real(1.0);
        real y = rndReal(rnd) * real(2.0) - real(1.0);
        real z = rndReal(rnd) * real(2.0) - real(1.0);
        
        real d = x*x+y*y+z*z;
        if(d > real(1)) {
            --i;
            continue;
        }
        
        p.pos.data[0] = x * diam + rndReal(rnd);
        p.pos.data[1] = y * diam + rndReal(rnd);
        p.pos.data[2] = z * diam + rndReal(rnd);

        real velVal = Math::sqrt(GRAV_CONST * centermass) * real(2.0);
        p.vel.data[0] = -p.pos.data[1] * velVal;
        p.vel.data[1] = p.pos.data[0] * velVal;
        p.vel.data[2] = p.pos.data[0] * velVal;
    }
}

static int simulate(particle *particles, int n, const real dt, const real theta, const AABB& maxAABB, const real maxAABBscaling)
{
	node *root = 0;

    #ifdef USE_BRUTE_FORCE_ALGO
        (void)theta;
        #ifdef USE_OPENMP
            #pragma omp parallel for schedule(dynamic, 10)
        #endif
        for(int i = 0; i < n; ++i) {
            particle& p = particles[i];
            for(int k = i+1; k < n; ++k) {
                particle& other = particles[k];
                addForce(p, other);
            }
        }
    #else
		root = addToTree(particles, n);
        #ifdef USE_OPENMP
            #pragma omp parallel for schedule(dynamic, 10)
        #endif
        for(int i = 0; i < n; ++i) {
            calcTreeForce(root, particles[i], theta);
        }
    #endif


    paint(particles, n, root);

    for(int i = 0; i < n; ++i) {
        integrate(particles[i], dt);
    }

    #ifdef USE_BRUTE_FORCE_ALGO
	#else
    destroyTree(root);
	#endif

    //remove particles
    for(int i = 0; i < n; ++i) {
        if(maxAABB.isInside(particles[i].pos, maxAABBscaling)) {
            continue;
        }
        if((n != 1) || (i != (n - 1))) {
            particles[i] = particles[n-1];
        }
        --n;
    }

    return n;
}

static void scaleScene(particle *particles, int n, real factor)
{
    for(int i = 0; i < n; ++i) {
        particle& p = particles[i];
        p.mass *= factor;
        p.mass *= factor;
        p.mass *= factor;
        p.pos = p.pos * factor;
        p.vel = p.vel * factor;
    }
}

int main(int argc, char **argv)
{
    int w = 800;
    int h = 800;
    int n = 20000;
    uint32_t seed = 1337;
    real dt = 2.0;
    real theta = 0.5;
    real maxAABBscaling = 8.0;
    real scaling = 0.1;

	(void)argc;
	(void)argv;
    //if(argc >= 9) {
    //    w = atoi(argv[1]);
    //    h = atoi(argv[2]);
    //    n = atoi(argv[3]);
    //    seed = atoi(argv[4]);
    //    dt = atof(argv[5]);
    //    theta = atof(argv[6]);
    //    maxAABBscaling = atof(argv[7]);
    //    scaling = atof(argv[8]);
    //}

    XSRandom rnd(seed);

    particle *particles = new particle[n];

    createScene(particles, n, rnd);
    scaleScene(particles, n, scaling);

    AABB aabb;
    for(int i = 0; i < n; ++i) {
        particle& p = particles[i];
        aabb.addVec(p.pos);
    }

    makeCube(aabb);

    initPaint(w, h, aabb);

    for(;;) {
        n = simulate(particles, n, dt, theta, aabb, maxAABBscaling);
    }

    delete [] particles;

    quitPaint();

    return 0;
}
