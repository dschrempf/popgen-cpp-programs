/**
 * @file   general_discrete_distributions.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Mar 18 11:11:35 2015
 * 
 * @brief Given K discrete events with different probabilities P[k],
 * produce a random value k consistent with its probability.
 * 
 * 
 */

#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main ()
{
    const gsl_rng_type * T;
    gsl_rng * R;

    size_t K = 4;
    const double P[4] = {0.1, 0.2, 0.3, 0.4};

    gsl_ran_discrete_t * G;

    size_t n;
    int N = 100000;
    int counter[4] = {0, 0, 0, 0};
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc (T);    
    gsl_rng_set (R, 1);

    G = gsl_ran_discrete_preproc (K, P);

    for (int i = 0; i < N; i++) {
        n = gsl_ran_discrete (R, G);
        counter[n]++;
    }

    for (int i = 0; i < 4; i++) {
        std::cout << (double) counter[i] /  (double) N << " ";
    }
    std::cout << std::endl;

    gsl_ran_discrete_free (G);
    gsl_rng_free (R);

    return 0;
}
