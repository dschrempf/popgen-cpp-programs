#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * @file   genetic_drift.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Feb 11 14:20:44 2015
 * 
 * @brief  Simulate genetic drift.
 *
 * Simulate genetic drift using the Wright-Fisher model; the GSL
 * library is used for random number generation.  Apart from that the
 * code is pretty straight forward to read.
 * 
 */


unsigned int
evolve_population (unsigned int k, unsigned int N, const gsl_rng * r)
{/** 
  * Evolve a population one Wright-Fisher generation.
  * 
  * @param k Number of individuals with allele one.
  * @param N Population size.
  * @param r GSL random number generator.
  * 
  * @return Number of individuals with allele one after one
  * Wright-Fisher generation.
  */
  double p;

  p = ((double) k) / ((double) N);
  
  return gsl_ran_binomial (r, p, N);
}

int
main (void)
{
  unsigned int N = 1000;        /* Total population size. */
  unsigned int k = 500;         /* Initial number of individuals with
                                   allele one. */

  unsigned int i = 0;

  const gsl_rng_type * T;
  gsl_rng * r;

  /* Create a generator chosen by the environment variable
     GSL_RNG_TYPE with seed GSL_RNG_SEED. */
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  std::cout << "Simulation of genetic drift." << std::endl;
  std::cout << "Starting values: k = " << k << ", N = " << N << std::endl;
  std::cout << k << std::endl;
  while ((k != 0) && (k != N)) {
    k = evolve_population (k, N, r);
    i++;
    std::cout << k << std::endl;
  }

  std::cout << "Total number of generations to fixation: " << i << std::endl;
  gsl_rng_free (r);

  return 0;
}
