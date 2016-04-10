#include "ctmc.h"

void set_q_no_diagonal(gsl_matrix * q, gsl_matrix * q_no_diagonal,
                  size_t ns) {
    size_t i, j;
    int above_diagonal;
    for (i = 0; i < ns; i++) {
        above_diagonal = 0;
        for (j = 0; j < ns-1; j++) {
            if (i == j) above_diagonal = 1;
            gsl_matrix_set(q_no_diagonal,
                           i, j,
                           gsl_matrix_get(q,i,j+above_diagonal));
        }
    }
}

CTMC::CTMC(gsl_matrix * q, size_t ns,
           bool log_path):
    q(q),
    ns(ns),
    log_out(std::cout)
{
    // General CTMC variables and parameters.
    time_now = 0;
    time_previous = 0;
    state_now = 0;
    state_previous = 0;
    q_no_diagonal = gsl_matrix_alloc (ns, ns-1);
    set_q_no_diagonal(q, q_no_diagonal, ns);

    // Auxiliary variables.
    rg = new RanGen();
    burn_in = BURN_IN;
    if (burn_in <= 100)
        out_warning("Low burn in value.");
    jump_counter = 0;

    // Analysis variables.
    invariant_distribution = new_zeroes<double>(ns);
    direct_hitting_times = gsl_matrix_alloc(ns, ns);
    times_direct_i = gsl_matrix_alloc(ns, ns);
    times_direct_n = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(direct_hitting_times);
    gsl_matrix_set_zero(times_direct_i);
    gsl_matrix_set_zero(times_direct_n);
    hitting_times = gsl_matrix_alloc(ns, ns);
    times_hitting_i = gsl_matrix_alloc(ns, ns);
    times_hitting_n = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(hitting_times);
    gsl_matrix_set_all(times_hitting_i, -1.0);
    gsl_matrix_set_zero(times_hitting_n);
    direct_number_jumps = gsl_matrix_alloc(ns, ns);
    jumps_direct_i = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(direct_number_jumps);
    gsl_matrix_set_zero(jumps_direct_i);

    // Logs.
    set_log_level();
    log_path = log_path;

    // Burn chain in.
    burn_it_in();
}

CTMC::~CTMC() {
    delete rg;
    gsl_matrix_free(q_no_diagonal);
    delete[] invariant_distribution;
    gsl_matrix_free(direct_hitting_times);
    gsl_matrix_free(times_direct_i);
    gsl_matrix_free(times_direct_n);
    gsl_matrix_free(hitting_times);
    gsl_matrix_free(times_hitting_i);
    gsl_matrix_free(times_hitting_n);
    gsl_matrix_free(direct_number_jumps);
    gsl_matrix_free(jumps_direct_i);
    if (log_path) {
        state_vector.clear();
        time_vector.clear();
    }
}


void CTMC::jump_maybe(double dt_max) {
    double holding_time =
        rg->pick_exponential(-1.0/gsl_matrix_get(q,state_now,state_now));
    if (holding_time <= dt_max) {
        time_previous = time_now;
        time_now += holding_time;
        jump();
        analyze_jump();
    }
    else time_now += dt_max;
}

void CTMC::jump() {
    state_previous = state_now;
    gsl_vector_view row = gsl_matrix_row(q_no_diagonal,state_now);
    state_now = rg->vector_weighted_pick(&row.vector,ns-1);
    // If the jump is not to the left, move one further up, because
    // the diagonal entry has been removed from the row,
    if (state_now >= state_previous) state_now++;
    if (log_path) log_push_back();
    jump_counter++;
    if (log_l >= 1) {
        if (jump_counter % LOG_INTERVAL == 0 ||
            log_l >= 3) {
            log_out << std::setiosflags(std::ios::fixed);
            log_out << "Jump ";
            log_out << std::setprecision(0) << std::setw(8);
            log_out << jump_counter;
            log_out << "; time ";
            log_out << std::setprecision(1) << std::setw(10);
            log_out << time_now;
            log_out << "; ";
            log_out << std::setprecision(0) << std::setw(4);
            log_out << state_previous;
            log_out << " to ";
            log_out << std::setprecision(0) << std::setw(4);
            log_out << state_now;
            log_out << "." << std::endl;
        }
    }
}

void CTMC::jump_silently() {
    state_previous = state_now;
    gsl_vector_view row = gsl_matrix_row(q_no_diagonal,state_now);
    state_now = rg->vector_weighted_pick(&row.vector,ns-1);
    if (state_now >= state_previous) state_now++;
}

void CTMC::burn_it_in() {
    unsigned int i;
    for (i = 0; i < burn_in; i++) jump_silently();
    if (log_l >= 1) {
        log_out << "Burn in of " << burn_in;
        log_out << " jumps finished." << std::endl;
    }
}

state CTMC::run(double dt) {
    double t_end = time_now + dt;
    while (time_now < t_end)
        jump_maybe(t_end - time_now);
    if (log_path) log_push_back();
    // The chain has come to an end :(.  Finish up the analysis.
    // Scale the invariant distribution and the matrizes.
    unsigned int i;
    for (i = 0; i < ns; i++) invariant_distribution[i] /= t_end;
    gsl_matrix_div_elements(direct_hitting_times, times_direct_n);
    gsl_matrix_div_elements(hitting_times, times_hitting_n);
    gsl_matrix_div_elements(direct_number_jumps, times_direct_n);
    return state_now;
}

void print_matrix(gsl_matrix * q, size_t s1, size_t s2,
                  std::ostream& out) {
    out << std::setiosflags(std::ios::fixed);
    out << std::setprecision(4);
    size_t i, j;
    for (i = 0; i < s1; i++) {
        for (j = 0; j < s2; j++) {
            out << std::setw(7) << gsl_matrix_get(q, i, j) << " ";
        }
        out << std::endl;
    }
}

void CTMC::print_info(std::ostream & out) {
    out << "------------------------------------------------------------" << std::endl;
    out << "Number of states: " << ns <<std::endl;
    out << "Current state: " << state_now << std::endl;
    out << "Current time: " << time_now << std::endl;
    out << "Transition rate matrix: " << std::endl;
    print_matrix(q, ns, ns, out);
    if (log_l >= 3) {
        out << "Transition rate matrix without diagonals: " << std::endl;
        print_matrix(q_no_diagonal, ns, ns-1, out);
    }
}

void CTMC::print_info_short(std::ostream& out) {
    out << "Current state: " << state_now << std::endl;
    out << "Current time: " << time_now << std::endl;
}

void CTMC::set_log_level(unsigned int level) {
    log_l = level;
}

void CTMC::log_push_back() {
    state_vector.push_back(state_now);
    time_vector.push_back(time_now);
}

void CTMC::print_path(std::ostream& out) {
    out << std::endl;
    out << std::setw(8) << "time";
    out << std::setw(8) << "state" << std::endl;;
    out << std::setiosflags(std::ios::fixed);
    std::vector<unsigned int>::iterator it1;
    std::vector<double>::iterator it2;
    for (it1 = state_vector.begin(), it2 = time_vector.begin();
         it1 != state_vector.end() || it2 != time_vector.end();
         it1++, it2++) {
        out << std::setprecision(3);
        out << std::setw(8) << *it2;
        out << std::setprecision(0);
        out << std::setw(8) << *it1 << std::endl;
    }
    out << std::endl;
}

void CTMC::analyze_jump() {
    double time_start;
    double time_delta;

    int jump_start;
    unsigned int jump_now = jump_counter;
    unsigned int jump_delta;
    
    unsigned int i = 0;
    
    // Increment invariant distribution.  It looks back one time
    // point.
    invariant_distribution[state_previous] +=
        (time_now - time_previous);
    // Set time_i.
    gsl_matrix_set_row_to_value(times_direct_i, state_now, time_now);
    gsl_matrix_set_row_to_value(jumps_direct_i, state_now, jump_now);
    // For the moving times we only reset the time if we have visited
    // state i inbetween.
    for (i = 0; i < ns; i++) {
        if (gsl_matrix_get(times_hitting_i, state_now, i) == -1.0)
            gsl_matrix_set(times_hitting_i, state_now, i, time_now);
    }
    
    // If initial time or jump is -1 then the initial state has not
    // been visited inbetween.
    for (i = 0; i < ns; i++) {
        if ((time_start =
             gsl_matrix_get(times_direct_i, i, state_now)) != -1.0) {
            time_delta = time_now - time_start;
            gsl_matrix_entry_add(direct_hitting_times,
                                 i, state_now, time_delta);
            gsl_matrix_entry_add(times_direct_n, i, state_now, 1.0);
            // Set the initial time to -1 to indicate that we used
            // this value.
            gsl_matrix_set(times_direct_i, i, state_now, -1.0);
        }
        if ((time_start =
             gsl_matrix_get(times_hitting_i, i, state_now)) != -1.0) {
            time_delta = time_now - time_start;
            gsl_matrix_entry_add(hitting_times,
                                 i, state_now, time_delta);
            gsl_matrix_entry_add(times_hitting_n, i, state_now, 1.0);
            gsl_matrix_set(times_hitting_i, i, state_now, -1.0);
        }
        if ((jump_start =
             gsl_matrix_get(jumps_direct_i, i, state_now)) != -1) {
            jump_delta = jump_now - jump_start;
            gsl_matrix_entry_add(direct_number_jumps,
                                 i, state_now, (double) jump_delta);
            // Same here for the initial jump.
            gsl_matrix_set(jumps_direct_i, i, state_now, -1);
        }
    }
}

void CTMC::print_direct_hitting_times(std::ostream& out) {
    out << "Average direct hitting times." << std::endl;
    gsl_matrix_print(direct_hitting_times, out);
}

void CTMC::print_hitting_times(std::ostream& out) {
    out << "Average moving times." << std::endl;
    gsl_matrix_print(hitting_times, out);
}

void CTMC::print_direct_number_jumps(std::ostream& out) {
    out << "Average number of jumps." << std::endl;
    gsl_matrix_print(direct_number_jumps, out);
}

void CTMC::print_invariant_distribution(std::ostream& out) {
    out << "The invariant distribution is:" << std::endl;
    unsigned int i;
    out << std::setiosflags(std::ios::fixed);
    out << std::setprecision(8);
    for (i = 0; i < ns; i++) {
        out << std::setw(10) << invariant_distribution[i] << std::endl;
    }
}
