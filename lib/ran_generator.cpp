#include "ran_generator.h"

RanGen::RanGen ()
{
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
}


RanGen::RanGen(unsigned long int s)
{
    RanGen();
    gsl_rng_set(r, s);
}
RanGen::~RanGen ()
{
    gsl_rng_free (r);
}


void RanGen::set_seed(unsigned long int s) {
    gsl_rng_set(r, s);
}

double RanGen::pick_exponential (double mean)
{
    return gsl_ran_exponential (r, mean);
}

double RanGen::pick_uniform()
{
    return gsl_rng_uniform (r);
}

int RanGen::vector_weighted_pick (gsl_vector * v, int l)
{
    double sum = 0.0;
    for (int i = 0; i < l; i++) sum += gsl_vector_get(v, i);

    double rn = pick_uniform ();
    double x = 0.0;
    for (int i = 0; i < l; i++) {
        x += gsl_vector_get(v,i);
        if (rn*sum < x) return i;
    }
    throw "Nothing has been picked.";
}
