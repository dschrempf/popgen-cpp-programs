/**
 * @file   ran_generator.h
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Tue Sep 22 17:19:29 2015
 *
 * @brief  A class to generate random variables.
 *
 *
 */

#ifndef RAN_GENERATOR_H
#define RAN_GENERATOR_H

#include <iomanip>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * The class RanGen is a wrapper around the GSL library.  It contains
 * useful functions to pick random variables out of different
 * distributions.
 * 
 */

class RanGen
{
 public:
    RanGen();
    ~RanGen();

    /** 
     * Simulate an exponentially distributed random variable.
     * 
     * @param mean the mean.
     * 
     * @return the picked value.
     */
    double pick_exponential (double mean);

    /** 
     * Simulate a uniformly distributed random variable between 0 and 1.
     * 
     * 
     * @return the picked value.
     */
    double pick_uniform ();

    /** 
     * Randomly pick an element out of a vector according to its value.
     *
     * For example, a for vector v = {0, 1, 1, 2}, v[0] will never be
     * chosen, v[3] will be chosen twice as often as v[1] or v[2]. 
     * 
     * @param v the vector.
     * @param l the length.
     * 
     * @return the picked index.
     */
    int vector_weighted_pick (gsl_vector * v, int l);

 private:
    const gsl_rng_type * T;
    gsl_rng * r;
};

#endif
