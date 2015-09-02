#include <iostream>
#include <gsl/gsl_rng.h>

/**
 * @file   ehrenfest_mcmc.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Feb 11 13:56:44 2015
 *
 * @brief Simulate gas particles in a divided box
 *
 * Please see Markov Chains --- Norris, exercise 1.7.2 (solution in
 * Feres, homework 6).
 * 
 * N gas particles move around in a box that is divided into two
 * symmetrical parts A (i particles) and B (N-i particles).  A hole
 * lets particles move from A to B and vice versa.  If we neglect two
 * particles passing at once, this model behaves like a Markov chain.
 * Here, we simulate this chain.
 *
 * The transition probabilites are as follows:
 * \image html ehrenfest_mcmc_transition_probabilities.png
 *
 * A simulation result:
 * \image html ehrenfest_mcmc_simulation.png
 * 
 */

unsigned int N = 100;           /**< Number of particles in the box. */
unsigned int i_start = 0;       /**< Starting state. */
unsigned int n_iter = 1000;     /**< Number of iterations to be simulated. */

// double get_p (unsigned int i, unsigned int N)
// {
//     return (double) (N-i) / (double) N;
// }

/**
 * Calculate the probability to jump from i to i-1.
 *
 * @param i Current state.
 * @param N Total number of particles in the box.
 *
 * @return Probability to jump from i to i-1.
 */
double get_q (unsigned int i, unsigned int N)
{
    return (double) i / (double) N;
}

/**
 * Let the chain jump from i to either i+1 or i-1.
 *
 * @param i Current state.
 * @param N Total number of particles in the box.
 * @param r GSL random number generator.
 *
 * @return New state.
 */
unsigned int jump (unsigned int i, unsigned int N, gsl_rng * r)
{
    // double p = get_p (i, N);
    double q = get_q (i, N);

    double s = gsl_rng_uniform (r);
    if (s < q) {
        return i-1;
    }
    else {
        return i+1;
    }
}

int main(void)
{
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    unsigned int states[n_iter];
    states[0] = i_start;
    std::cout << states[0] << std::endl;

    for (unsigned int j = 1; j < n_iter; j++) {
        states[j] = jump (states[j-1], N, r);
        std::cout << states[j] << std::endl;
    }

    gsl_rng_free (r);
    return 0;
}
