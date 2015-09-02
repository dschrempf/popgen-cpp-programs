#include <iostream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * @file   continuous_markov_chain_norris_ex_2_3_2.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Sep  2 12:21:18 2015
 * 
 * @brief  Runs a continuous Markov chain.
 * 
 * Please see Markov Chains --- Norris, exercise 2.3.2 (solution in
 * Feres, homework 8).
 * 
 * The problem is the following: we draw independent exponential
 * random variables \f$ T_1, T_2, \ldots \f$ of parameter \f$ \lambda
 * \f$.  We do this \f$ N \f$ times, where \f$ N \f$ is an independent
 * geometric random variable of parameter \f$ \beta \f$.
 *
 * The hypothesis is: the stopping time \f$ T = T_1 + T_2 + \ldots +
 * T_N \f$ after these draws has exponential distribution of parameter
 * \f$ \lambda \beta \f$.
 *
 * This can be proven analytically (see the notes of Feres) and can
 * also be made visible by simulations.  Here, we see that for \f$
 * \lambda=1.0 \f$, the mean stopping time is just \f$ 1/\beta \f$.
 * This is exactly the mean of the exponential distribution with
 * parameter \f$ \beta \f$.
 * \image html continuous_markov_chain_norris_ex_2_3_2.png
 * 
 */

double lambda = 1.0;            /**< Parameter of exponential
                                   distribution. */
int n_iter    = 1000;           /**< Number of iterations. */
int n_points  = 100;            /**< Number of grid points for plot of
                                   exponential distribution of
                                   parameter beta*lambda. */

class Chain
{
public:
    Chain () {t=0; s=0;};
    Chain (double t_init, int s_init) {t=t_init; s=s_init;};
    double t;                   /**< time */
    int s;                      /**< state */
};

bool jump (Chain * chain, double beta, gsl_rng * r)
{
    double holding_time = gsl_ran_exponential(r, 1.0/lambda);
    double d = gsl_rng_uniform (r);
    chain->s+=1;
    chain->t+=holding_time;
    if (d < beta) return false;
    return true;
}

int main(void)
{
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    double sigma[n_points];     /**< Mean state upon absorption. */
    double tau[n_points];       /**< Mean time upon absorption. */

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    std::cout << std::setw(7) << "beta";
    std::cout << std::setw(7) << "sigma";
    std::cout << std::setw(7) << "tau" << std::endl;
    for (int b = 1; b <= n_points; b++) {
        Chain chain;
        double beta = (double) b / (double) n_points;
        for (int i = 0; i < n_iter; i++) {
            while (jump (&chain, beta, r));      
        }
        sigma[b] = chain.s / (double) n_iter;
        tau[b] = (double) chain.t / (double) n_iter;
        std::cout << std::setw(7) << beta;
        std::cout << std::setw(7) << sigma[b];
        std::cout << std::setw(7) << tau[b] << std::endl;
    }

    return 0;
}
