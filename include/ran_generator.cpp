#include "ran_generator.h"

RanGen::RanGen ()
{
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
}

RanGen::~RanGen ()
{
    gsl_rng_free (r);
}

double RanGen::get_ran_exponential (double mean)
{
    return gsl_ran_exponential (r, mean);
}

double RanGen::get_uniform()
{
    return gsl_rng_uniform (r);
}

int RanGen::ran_pick (double * v, int l)
{
    double sum = 0.0;
    for (int i = 0; i < l; i++) sum += v[i];

    double rn = gsl_rng_uniform (r);
    double x = 0.0;
    for (int i = 0; i < l; i++) {
        x += v[i];
        if (rn*sum < x) return i;
    }
    throw "Nothing has been picked.";
}
