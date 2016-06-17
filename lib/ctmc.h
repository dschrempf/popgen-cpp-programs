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

typedef unsigned int state;

unsigned int BURN_IN = 10000;
unsigned int LOG_INTERVAL = 1e5;

class CTMC {
 public:
    /** Initialize the chain.
     *
     *  @param q the transition rate matrix Q.
     *  @param ns the number of states.
     *  @param log_path set to true to log the full path.
     */
    CTMC(gsl_matrix * q, size_t ns,
         bool log_path=false);

    ~CTMC();

    /**
     * Let the chain jump once if the time to the next jump does not
     * exceed the given value.
     *
     * Wait some time t that is exponentially distributed with 1/q_ii.
     * If t < dt_max, let the chain jumps once.
     *
     */
    void jump_maybe(double dt_max);

    /**
     * Let the chain jump one.
     *
     * The time is not changed, the new state is determined according
     * to the distribution of the q_ij.
     *
     */
    void jump();

    /**
     * Same as jump() but not logs, etc.
     *
     */
    void jump_silently();

    // TODO: There might be a bug because i do not honor the expected
    // times to the next jump.
    /**
     * Evolve the chain for the number of jumps specified in burn_in.
     *
     */
    void burn_it_in();

    /**
     * Let the chain run for the given time.
     *
     * @param dt the time to run.
     *
     * @return the state the chain ends up in.
     */
    state run(double dt);

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
     *nnn
     *
     * @param level log level.
     */
    void set_log_level(unsigned int level=1);

    /**
     * Log current state of Markov chain.
     *
     */
    void log_push_back();

    /**
     * Print the states and times.
     *
     * @param out output stream.
     */
    void print_path(std::ostream & out);

    const std::vector<unsigned int> get_state_vector() const
    {
    return state_vector;
    }

    const std::vector<double> get_time_vector() const
    {
    return time_vector;
    }

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
    void analyze_jump();

    void print_direct_hitting_times(std::ostream & out);

    void print_hitting_times(std::ostream & out);

    void print_direct_number_jumps(std::ostream & out);

    void print_invariant_distribution(std::ostream & out);

    double get_entry_invariant_distribution(unsigned int i);
    
    /// Random number generator.
    RanGen * rg;
    
 private:
    /// Time of the Markov chain.
    double time_now;
    /// Time of previous jump of the Markov chain.
    double time_previous;
    /// Current state of the chain.
    state state_now;
    /// Previous state of the chain.
    state state_previous;
    /// Transition rate matrix Q.
    gsl_matrix * q;
    /// Number of states.
    size_t ns;
    /// Transition rate matrix Q without diagonal elements.
    /// Internally used to determine the state to jump to.
    gsl_matrix * q_no_diagonal;
    /// Counter of jumps.
    unsigned int jump_counter;
    /// Number of burn in jumps.
    unsigned int burn_in;
    /// Debug level (0 to 3).
    unsigned int log_l;
    /// The output stream to write logs to.
    std::ostream& log_out;
    /// Log states and times?
    bool log_path;
    /// Logged states.
    std::vector<unsigned int> state_vector;
    /// Logged times.
    std::vector<double> time_vector;
    /// The average direct hitting times.  The entry (i,j) is the
    /// average direct hitting time from state i to state j.  I.e.,
    /// the time to move from state i to state j without moving back
    /// to i inbetween.
    gsl_matrix * direct_hitting_times;
    /// Initial times; direct hitting times analysis.
    gsl_matrix * times_direct_i;
    /// Number of hits matrix; direct hitting times analysis.
    gsl_matrix * times_direct_n;
    /// Hitting time.  I am interested in the time that the chain
    /// needs from state i to state j.  This is the true hitting time
    /// and differs from the direct hitting time.  It is the time FROM
    /// - the first visit to i after a visit to j TO
    /// - the next visit to j.
    gsl_matrix * hitting_times;
    gsl_matrix * times_hitting_i;
    gsl_matrix * times_hitting_n;
    /// The average number of jumps.  The entry (i,j) is the average
    /// number of jumps from state i to state j.
    gsl_matrix * direct_number_jumps;
    /// Initial jump numbers; number of jumps analysis.
    gsl_matrix * jumps_direct_i;
    /// The invariant distirbution.
    double * invariant_distribution;
};

#endif
