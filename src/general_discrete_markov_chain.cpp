/**
 * @file   general_discrete_markov_chain.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Mar 18 15:51:19 2015
 *
 * @brief Simulate a general discrete Markov chain with a given transition
 * probability matrix P.
 *
 * Assume that P has a stationary distribution pi.  Iterate the chain
 * and print the proportion of time spent in each state.  This should
 * match the stationary distribution.
 *
 * @todo Object oriented programming.
 * 
 */

#include <iostream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//////////////////////////////////////////////////////////////////////
// User interface.
size_t n_states       = 4;      /**< Number of states. */
unsigned int s_init   = 0;      /**< Initial state. */
/// Transition probability matrix.
double trans_matrix[] = { 0.1,  0.4,  0.4,  0.1,
                          0.25, 0.25, 0.25, 0.25,
                          0.01, 0.01, 0.97, 0.01,
                          0.33, 0.33, 0.33, 0.01};
//////////////////////////////////////////////////////////////////////
// Markov chain settings.
unsigned int n_iter = 10000000;

//////////////////////////////////////////////////////////////////////
/**
 * Perform a single jump of the discrete Markov chain.
 * 
 * @param G A GSL discrete random number lookup table for the current
 * state.
 * @param R A random number generator.
 * 
 * @return New state.
 */
int jump(const gsl_ran_discrete_t * G, const gsl_rng * R)
{
    return gsl_ran_discrete (R, G);
}

void pp_matrix (double * M, size_t r, size_t c) {
    for (unsigned int i = 0; i < r; i++) {
        for (unsigned int j = 0; j < c; j++) {
            std::cout
                << std::setprecision(2)
                << std::fixed
                << std::setw(4)
                << M[i*c + j] << " ";
        }
        std::cout << std::endl;
    }
}

int main ()
{
    unsigned int counter_l[n_states];
    unsigned int s = s_init;
    const gsl_rng_type * T;
    gsl_rng * R;

    gsl_ran_discrete_t * G[n_states];
    const double * P[n_states];

    // Setup counter.
    for (unsigned int i = 0; i < n_states; i++)
        counter_l[i] = 0;

    // Setup random number generator.
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc(T);
    gsl_rng_set (R, 1);

    // Setup discrete random number generator.
    for (unsigned int i = 0; i < n_states; i++) {
        P[i] = &trans_matrix[i*n_states];
        G[i] = gsl_ran_discrete_preproc (n_states, P[i]);
    }

    // Output information.
    std::cout << "The transition probability matrix P is:" << std::endl;
    pp_matrix(trans_matrix, n_states, n_states);
    std::cout << std::endl;
    std::cout << "The starting state is: " << s_init << std::endl;

    // Run the chain.
    for (unsigned int i = 0; i < n_iter; i++) {
        s = jump(G[s], R);
        counter_l[s]++;
    }

    // Summarize results.
    std::cout << n_iter << " iterations of the Markov chain have been performed."
              << std::endl;
    std::cout << "The proportion of time t spent in each state s is:"
              << std::endl << std::endl;
    std::cout << "s\tt" << std::endl;
    for (unsigned int i = 0; i < n_states; i++) {
        double prop = (double) counter_l[i] / (double) n_iter;
        std::cout
            << std::setprecision(3)
            << i <<"\t" << prop << std::endl;
    }
    std::cout << std::endl;
    std::cout << "This should match the stationary distribution of P."
              << std::endl;

    // Free memory.
    for (unsigned int i = 0; i < n_states; i++)
        gsl_ran_discrete_free (G[i]);
    gsl_rng_free (R);
    return 0;
}
