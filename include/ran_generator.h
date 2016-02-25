/**
 * @file   ran_generator.h
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Tue Sep 22 17:19:29 2015
 *
 * @brief  A class to generate random variables.
 *
 *
 */

#include <iostream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class RanGen
{
 public:
    RanGen();
    ~RanGen();
    double get_ran_exponential (double mean);
    double get_uniform ();
    int ran_pick (double * v, int l);
 private:
    const gsl_rng_type * T;
    gsl_rng * r;
};
