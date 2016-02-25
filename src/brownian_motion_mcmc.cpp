/**
 * @file   brownian_motion_mcmc.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Thu Feb 12 13:57:43 2015
 *
 * @brief  Simulate standard Brownian motion (Wiener process).
 *
 * <a href="/home/dominik/Shared/university/E-Books/Mathematics/Topics in Applied Mathematics - Feres - Lecture Notes/dhigham.pdf">(Higham 2001)</a>
 * 
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(void)
{
    const gsl_rng_type * rng_type;
    gsl_rng * r;
    gsl_rng_env_setup ();
    rng_type = gsl_rng_default;
    r = gsl_rng_alloc (rng_type);

    double T  = 1;              // Simulation end time.
    int N     = 500;            // Number of iterations.
    double dt = T / (double) N; // Discrete time step size.
    double dW;                  // Discrete step size after dt.
    double W[N+1];              // Vector of particle positions.

    for (int i = 0; i < N; i++) {
        W[i] = 0;
    }

    double sqrtdt = sqrt(dt);
    T = 0;
    std::cout << std::fixed;
    std::cout << std::setw(12) << std::setprecision(4);
    std::cout << T << "\t";
    std::cout << std::setw(12) << std::setprecision(8);
    std::cout << W[0] << std::endl;
    for (int i = 1; i <= N; i++) {
        dW = sqrtdt * gsl_ran_gaussian (r, 1.0);
        W[i] = W[i-1] + dW;
        T += dt;
        std::cout << std::setw(12) << std::setprecision(4);
        std::cout << T << "\t";                                             
        std::cout << std::setw(12) << std::setprecision(8);
        std::cout << W[i] << std::endl;
    }

    gsl_rng_free (r);
    return 0;
}
