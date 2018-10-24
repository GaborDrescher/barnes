#ifndef XOR_SHIFT_RANDOM_INCLUDE
#define XOR_SHIFT_RANDOM_INCLUDE

#include "config.h"

// this uses xorshift1024, see https://en.wikipedia.org/wiki/Xorshift

class XSRandom
{
    private:
    uint64_t s[16];
    int p;

    public:
    XSRandom(uint64_t se = 42) {
        p = 0;
        seed(se);
    }
    void seed(uint64_t x)
    {
        //seed with xorshift64star
        for(int i = 0; i < 16; ++i) {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            s[i] = x * 2685821657736338717ull;
        }
    }

    uint64_t next()
    { 
        uint64_t s0 = s[p];
        uint64_t s1 = s[p = ( p + 1 ) & 15];
        s1 ^= s1 << 31;
        s1 ^= s1 >> 11;
        s0 ^= s0 >> 30;
        return (s[p] = s0 ^ s1) * 1181783497276652981ull;
    }

    uint64_t max()
    {
        return UINT64_MAX;
    }
};

#endif
