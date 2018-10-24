#ifndef BHT_MATH_INCLUDE
#define BHT_MATH_INCLUDE

// real type
#define real float
#define REAL_MAX ( 1E38f)
#define REAL_MIN (-1E38f)

// physical constants
#define GRAV_CONST real(6.673848E-11)
#define NEAR_EPS real(0.1)

class Math
{
    public:
    template<typename T>
    static T min(T a, T b)
    {
        return a < b ? a : b;
    }

    template<typename T>
    static T max(T a, T b)
    {
        return a > b ? a : b;
    }

    static real sqrt(real val)
    {
        return __builtin_sqrt(val);
    }
};

#endif

