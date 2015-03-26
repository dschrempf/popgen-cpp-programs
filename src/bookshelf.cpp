/**
 * @file   bookshelf.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Thu Mar 12 17:49:27 2015
 *
 * @brief  Bookshelf Markov chain.
 *
 * See Markov Chains Norris Ex. 1.8.4
 *
 */

#include <iostream>
#include <iomanip>

#include <gsl/gsl_rng.h>

// User interface.
int n_books = 7;              /**< Number of books in the shelf. */

// Markov chain settings.
int n_burn = 10000;
int n_iter = 1000000;

void jump (int * state_vec, const gsl_rng * r) {
    int taken = (int) gsl_rng_uniform_int (r, n_books);
    int tmp = state_vec[0];
    state_vec[0] = state_vec[taken];
    state_vec[taken] = tmp;
}

bool check (int * state_vec) {
    for (int i = 0; i < n_books; i++) {
        if (state_vec[i] != i) return false;
    }
    return true;
}

unsigned long int fac (unsigned long int n) {
    if (n == 1) return 1;
    return n * fac(n-1);
}

int main ()
{
    // Setup random number generator.
    const gsl_rng_type * t;
    gsl_rng * r;

    gsl_rng_env_setup();
    t = gsl_rng_default;
    r = gsl_rng_alloc (t);
    gsl_rng_set (r, 1);

    int state_vec[n_books];
    for (int i = 0; i < n_books; i++) {
        state_vec[i] = i;
    }


    // Burn in phase.
    for (int i = 0; i < n_burn; i++) {
        jump(state_vec, r);
    }
    
    int counter = 0;
    // MCMC steps.
    for (int i = 0; i < n_iter; i++) {
        jump(state_vec, r);
        if (check(state_vec)) counter++;
    }

    // for (int i = 0; i < n_books; i++) {
    //     std::cout << state_vec[i] << std::endl;
    // }
    double rel_occ = (double) counter / (double) n_iter;
    double n_states = 1.0 / rel_occ;
    std::cout << "The chain returned " << counter << " times to its initial state." << std::endl;
    std::cout << "The relative occurrence is " << rel_occ << "." << std::endl;
    std::cout << "This corresponds to an empirical number of " << n_states << " states." << std::endl;
    std::cout << "The actual number of states is " << fac(n_books) << "."<< std::endl;

    gsl_rng_free (r);
    return 0;
}
