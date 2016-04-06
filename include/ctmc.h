/**
 * @file   ctmc.h
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Mon Apr  4 09:51:24 2016
 *
 * @brief  A class to run a continuous-time Markov chain (CTMC).
 *
 *
 */

#ifndef CTMC_H
#define CTMC_H

#include <iostream>
#include <gsl/gsl_matrix.h>
#include <vector>
#include "ran_generator.h"
#include "tools.h"

class CTMCPath;

typedef unsigned int state;

unsigned int BURN_IN = 1000;
unsigned int LOG_INTERVAL = 1e5;

class CTMC {
 public:
    /** Initialize the chain.
     *
     *  @param q_matrix the transition rate matrix Q.
     *  @param n_states the number of states.
     *  @param log_out the log output stream.
     *  @param do_log_path set to false to inhibit logging of path.
     */
    CTMC(gsl_matrix * q_matrix, size_t n_states,
         std::ostream& log_out, bool do_log_path=true);

    /**
     *  Initialize the chain with a specific time and state.
     *
     *  @param q_matrix the transition rate matrix Q.
     *  @param n_states the number of states.
     *  @param t_init the initial time.
     *  @param s_init the inital state.
     *  @param log_out the log output stream.
     *  @param do_log_path set to false to inhibit logging.
     */
    CTMC(gsl_matrix * q_matrix, size_t n_states,
         double t_init, state s_init,
         std::ostream& log_out, bool do_log_path=true);

    ~CTMC();

    /**
     * Let the chain jump once if the time to the next jump does not
     * exceed the given value.
     *
     * Wait some time t that is exponentially distributed with 1/q_ii.
     * If t < dt_max, let the chain jumps once.
     *
     * @return the new state the chain ends up in.
     */
    state jump_maybe(double dt_max);

    /**
     * Let the chain jump one.
     *
     * The time is not changed, the new state is determined according
     * to the distribution of the q_ij.
     *
     *
     * @return
     */
    state jump();

    /**
     * Let the chain evolve for some time.
     *
     * @param dt the time to evolve.
     *
     * @return the state the chain ends up in.
     */
    state evolve(double dt);

    /**
     * Print info about the CTMC.
     *
     * @param out output stream.
     */
    void print_info(std::ostream & out);

    /**
     * Print short info about the CTMC.
     *
     * @param out output stream.
     */
    void print_info_short(std::ostream & out);

    /**
     * Set the log level from 0 (silent) to 3 (debug).  Default is 1.
     *
     *
     * @param level log level.
     */
    void set_log_level(unsigned int level=1);

    /**
     * Log the path of the CTMC (log states and time in time_vector
     * and state_vector).
     *
     */
    void activate_log_path();

    /**
     * Log current state of Markov chain.
     *
     */
    void log_push_back();

    /**
     * Pretty print the states and times.
     *
     * @param out output stream.
     */
    void pretty_print_path(std::ostream & out);

    const std::vector<unsigned int> get_state_vector() const
    {
    return state_vector;
    }

    const std::vector<double> get_time_vector() const
    {
    return time_vector;
    }

    /// The log class of the continuous time Markov chain.
    CTMCPath * ctmc_path;

 private:
    /// Time of the Markov chain.
    double t;
    /// Current state of the chain.
    state s;
    /// Number of states.
    size_t ns;
    /// Transition rate matrix Q.
    gsl_matrix * q;
    /// Transition rate matrix Q without diagonal elements.
    /// Internally used to determine jumps.
    gsl_matrix * q_no_diagonal;
    /// Random number generator.
    RanGen * rg;
    /// Counter of jumps.
    unsigned int counter;
    /// Debug level (0 to 3).
    unsigned int logl;
    /// Reference to the log output stream; set in the constructor
    /// (defaults to cout).
    std::ostream& log_out;
    /// Log states and times?
    bool log_path;
    /// Logged states.
    std::vector<unsigned int> state_vector;
    /// Logged times.
    std::vector<double> time_vector;
};

class CTMCPath {
 public:
    /**
     * Enter states and times and the number of states.
     *
     * @param p_states_vector pointer to states vector.
     * @param p_times_vector pointer to times vector.
     * @param n_states number of states.
     * @param log_level the log level (0 - silent, 3 - debug).
     * @param log_out the log output stream.
     * @param burn_in_init number of burn in jumps.
     */
    CTMCPath(std::vector<unsigned int> * p_states_vector,
             std::vector<double> * p_times_vector,
             unsigned int n_states,
             unsigned int log_level,
             std::ostream& log_out,
             unsigned int burn_in_init=BURN_IN);
    /**
     * Constructor that gets information out of a given CTMC.
     *
     * @param chain a chain.
     */
    CTMCPath(CTMC * chain, unsigned int burn_in_init=BURN_IN);

    ~CTMCPath();

    /**
     * Traverses logs and analyzes them.
     * - Determins the average hitting time from state a to state b.
     * - Determins the average number of jumps from state a to state b.
     *
     * @param out_average_hitting_time OUT.
     * @param out_average_n_jumps OUT.
     * @param burn_in number of jumps to wait before checking times.
     *
     * @return average hitting time.
     */

    void analyze(unsigned int burn_in=BURN_IN);

    void print_hitting_times(std::ostream & out);

    void print_n_jumps(std::ostream & out);

    void print_invariant_distribution(std::ostream & out);

 private:
    /// Pointer to the states vector.
    std::vector<unsigned int> * p_states;
    /// Pointer to the times vector.
    std::vector<double> * p_times;
    /// Number of states.
    size_t ns;
    /// Log length.
    size_t log_len;
    /// Log level;
    unsigned int logl;
    /// Number of burn in jumps.
    unsigned int burn_in;
    /// Reference to the log output stream; set in the constructor
    /// (default: cout).
    std::ostream& log_out;
    /// The average hitting times.  The entry (i,j) is the average
    /// hitting time from state i to state j.
    gsl_matrix * hitting_times;
    /// The average number of jumps.  The entry (i,j) is the average
    /// number of jumps from state i to state j.
    gsl_matrix * n_jumps;
    /// The invariant distirbution.
    double * invariant_distribution;
};

#endif
