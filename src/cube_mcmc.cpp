#include <iostream>
#include <gsl/gsl_rng.h>

/**
 * @file   cube_mcmc.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Feb 11 18:33:21 2015
 * 
 * @brief Markov Chains - Norris; exercise 1.7.3.; simulation of a
 * Markov chain that moves around the eight vertices of a cube.
 * 
 * 
 * Solution in Feres, homework 6.
 *
 * Simulate a Markov chain that moves around the eight vertices of a
 * cube.  Each step, the chain is equally likely to move to each of
 * the three adjacent vertices, independently of its past move.  Let i
 * be the initial vertex occupied by the chain, o the vertex opposite
 * i:
 *
 * \image html cube_mcmc_states.png
 *
 * The following quantities are calculated:
 *
 * 1. the expected number of steps until the particle returns to i;
 * 2. the expected number of visits to o until the first return to i;
 * 3. the expected number of steps until the first visit to o.
 *
 */

typedef int state_t;
typedef int move_t;

state_t s_init  = 0;            /**< Starting state. */
state_t s_opposite = 6;         /**< State opposite from the starting state. */
int n_turns = 100000;           /**< Number of turns to average over. */
int n_directions = 3;           /**< Number of direction to move to. */

move_t get_next_move (const gsl_rng * r)
{
    return (move_t) gsl_rng_uniform_int (r, n_directions);
}

state_t get_next_state (state_t s, move_t m)
{
    if (s == 0) {
        if (m == 0) return (state_t) 1;
        if (m == 1) return (state_t) 3;
        if (m == 2) return (state_t) 4;
    }
    if (s == 1) {
        if (m == 0) return (state_t) 2;
        if (m == 1) return (state_t) 0;
        if (m == 2) return (state_t) 5;
    }
    if (s == 2) {
        if (m == 0) return (state_t) 3;
        if (m == 1) return (state_t) 1;
        if (m == 2) return (state_t) 6;
    }
    if (s == 3) {
        if (m == 0) return (state_t) 0;
        if (m == 1) return (state_t) 2;
        if (m == 2) return (state_t) 7;
    }
    if (s == 4) {
        if (m == 0) return (state_t) 5;
        if (m == 1) return (state_t) 7;
        if (m == 2) return (state_t) 0;
    }
    if (s == 5) {
        if (m == 0) return (state_t) 6;
        if (m == 1) return (state_t) 4;
        if (m == 2) return (state_t) 1;
    }
    if (s == 6) {
        if (m == 0) return (state_t) 7;
        if (m == 1) return (state_t) 5;
        if (m == 2) return (state_t) 2;
    }
    if (s == 7) {
        if (m == 0) return (state_t) 4;
        if (m == 1) return (state_t) 6;
        if (m == 2) return (state_t) 3;
    }
    throw "Invalid state.";
}

state_t jump (state_t s, const gsl_rng * r)
{
    move_t m       = get_next_move (r);
    state_t next_s = get_next_state (s, m);
    return next_s;
}

double expected_steps_return (const gsl_rng * r)
{
    state_t s;
    int counter;
    int total_counts = 0;
    double exp;
    
    for (int i = 0; i < n_turns; i++) {
        s = jump (s_init, r);
        counter = 1;
        while (s != s_init) {
            s = jump (s, r);
            counter++;
        }
        total_counts += counter;
    }
    exp = (double) total_counts / (double) n_turns;
    return exp;
}

double expected_visits_o (const gsl_rng * r)
{
    state_t s;
    int visits = 0;
    double exp;
    
    for (int i = 0; i < n_turns; i++) {
        s = jump (s_init, r);
        if (s == s_opposite) visits++;
        while (s != s_init) {
            s = jump (s, r);
            if (s == s_opposite) visits++;
        }
    }
    exp = (double) visits / (double) n_turns;
    return exp;
}

double expected_steps_to_o (const gsl_rng *r)
{
    state_t s;
    int steps = 0;
    int total_steps = 0;
    double exp;
    
    for (int i = 0; i < n_turns; i++) {
        s = jump (s_init, r);
        steps = 1;
        while (s != s_opposite) {
            s = jump (s, r);
            steps++;
        }
        total_steps += steps;
    }
    exp = (double) total_steps / (double) n_turns;
    return exp;    
}

int main(void)
{
    const gsl_rng_type * T;
    gsl_rng * r; 

    double exp_return;
    double exp_visits;
    double exp_steps_o;
    
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    exp_return = expected_steps_return (r);
    exp_visits = expected_visits_o (r);
    exp_steps_o = expected_steps_to_o (r);
    
    std::cout << "Expected number of jumps before returing to the inital state: ";
    std::cout << exp_return << std::endl;
    std::cout << "Expected number of visits of the opposite state: ";
    std::cout << exp_visits << std::endl;
    std::cout << "Expected number of steps to the opposite state: ";
    std::cout << exp_steps_o << std::endl;
    gsl_rng_free (r);
    return 0;
}
