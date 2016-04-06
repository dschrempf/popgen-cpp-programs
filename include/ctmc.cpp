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

CTMC::CTMC(gsl_matrix * q_matrix, size_t n_states,
           std::ostream& log_out, bool do_log_path):
    log_out(log_out)
{
    t = 0;
    s = 0;
    ns = n_states;
    q = q_matrix;
    // q_no_diagonal is used to determine the state to jump to.
    // Allocate a matrix and set it accordingly.
    q_no_diagonal = gsl_matrix_alloc (ns, ns-1);
    set_q_no_diagonal(q, q_no_diagonal, ns);
    rg = new RanGen();
    counter = 0;
    set_log_level();
    // Initialize log_path to false and avtivate it if specified.
    log_path = false;
    if (do_log_path) {
        activate_log_path();
    }
}

CTMC::CTMC(gsl_matrix * q_matrix, size_t n_states,
           double t_init, state s_init,
           std::ostream& log_out, bool do_log_path):
    log_out(log_out)
{
    CTMC(q_matrix, n_states, log_out, do_log_path);
    t = t_init;
    s = s_init;
}

CTMC::~CTMC() {
    delete rg;
    gsl_matrix_free(q_no_diagonal);
    if (log_path) {
        state_vector.clear();
        time_vector.clear();
        delete ctmc_path;
    }
}


state CTMC::jump_maybe(double dt_max) {
    double holding_time = rg->pick_exponential(-1.0/gsl_matrix_get(q,s,s));

    if (holding_time <= dt_max) {
        s = jump();
        t = t + holding_time;
    }
    else {
        t = t + dt_max;
    }

    return s;
}

state CTMC::jump() {
    state new_state;
    gsl_vector_view row = gsl_matrix_row(q_no_diagonal,s);
    new_state = rg->vector_weighted_pick(&row.vector,ns-1);
    // If the jump is not to the left, move one further up, because
    // the diagonal entry has been removed from the row,
    if (new_state >= s) {
        new_state++;
    }
    if (log_path) log_push_back();
    counter++;
    if (logl >= 1) {
        if (counter % LOG_INTERVAL == 0 ||
            logl >= 3) {
            log_out << std::setiosflags(std::ios::fixed);
            log_out << "Jump ";
            log_out << std::setprecision(0) << std::setw(8);
            log_out << counter;
            log_out << "; time ";
            log_out << std::setprecision(1) << std::setw(7);
            log_out << t;
            log_out << "; ";
            log_out << std::setprecision(0) << std::setw(4);
            log_out << s;
            log_out << " to ";
            log_out << std::setprecision(0) << std::setw(4);
            log_out << new_state;
            log_out << "." << std::endl;
        }
    }
    s = new_state;
    return s;
}

state CTMC::evolve(double dt) {
    double t_end = t + dt;
    while (t < t_end)
        jump_maybe(t_end - t);
    if (log_path) log_push_back();
    return s;
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
    out << "Current state: " << s << std::endl;
    out << "Current time: " << t << std::endl;
    out << "Transition rate matrix: " << std::endl;
    print_matrix(q, ns, ns, out);
    if (logl >= 3) {
        out << "Transition rate matrix without diagonals: " << std::endl;
        print_matrix(q_no_diagonal, ns, ns-1, out);
    }
}

void CTMC::print_info_short(std::ostream& out) {
    out << "Current state: " << s << std::endl;
    out << "Current time: " << t << std::endl;
}

void CTMC::set_log_level(unsigned int level) {
    logl = level;
}

void CTMC::activate_log_path() {
    if (log_path) out_warning("Logs have already been activated.");
    else {
        log_path = true;
        ctmc_path = new CTMCPath(&state_vector,
                                 &time_vector, ns,
                                 logl, log_out);
    }
}

void CTMC::log_push_back() {
    state_vector.push_back(s);
    time_vector.push_back(t);
}

void CTMC::pretty_print_path(std::ostream& out) {
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

CTMCPath::CTMCPath(std::vector<unsigned int> * p_states_vector,
                   std::vector<double> * p_times_vector,
                   unsigned int n_states,
                   unsigned int log_level,
                   std::ostream& log_out,
                   unsigned int burn_in_init):
    log_out(log_out)
{
    unsigned int size = p_states_vector->size();
    if (size != p_times_vector->size())
        throw "Vectors do not have equal size.";
    burn_in = burn_in_init;
    log_len = size;
    p_states = p_states_vector;
    p_times = p_times_vector;
    ns = n_states;
    logl = log_level;
    invariant_distribution = new_zeroes<double>(ns);
    // Initialize hitting times and number of jump matrizes.
    hitting_times = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(hitting_times);
    n_jumps = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(n_jumps);
}

CTMCPath::~CTMCPath() {
    delete[] invariant_distribution;
    gsl_matrix_free(hitting_times);
    gsl_matrix_free(n_jumps);
}

void CTMCPath::analyze(unsigned int burn_in) {
    // Initial times.
    gsl_matrix * times_i = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(times_i);
    // Number of hits matrix.
    gsl_matrix * times_n = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(times_n);
    double time_start;
    double time_now;

    // Initial jump numbers.
    gsl_matrix * jumps_i = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(jumps_i);
    // Number of hits matrix (using both times_n and jumps_n is
    // actually redundant, but easier to read).
    gsl_matrix * jumps_n = gsl_matrix_alloc(ns, ns);
    gsl_matrix_set_zero(jumps_n);
    unsigned int jump_start;
    unsigned int jump_now;

    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int c = 0;
    unsigned int state_previous;
    unsigned int state_now;
    
    // Error checking.
    if (burn_in >= p_states->size())
        throw "Burn in larger than length of chain.";
    if (burn_in <= 100)
        out_warning("Low burn in value.");
    if (burn_in <= 0)
        throw "Logging only possible after first jump.";

    // Traverse through CTMC path from burn_in to the end.
    std::vector<unsigned int>::iterator it;
    for (i = burn_in; i < p_states->size(); i++) {
        // Set the current state.
        state_previous = p_states->at(i-1);
        state_now = p_states->at(i);
        time_now = p_times->at(i);
        jump_now = i;
        // Increment invariant distribution.  It looks back one time
        // point.
        invariant_distribution[state_previous] +=
            (p_times->at(i) - p_times->at(i-1));
        // Set time_i.
        gsl_matrix_set_row_to_value(times_i, state_now, time_now);
        gsl_matrix_set_row_to_value(jumps_i, state_now, i);
        
        for (j = 0; j < ns; j++) {
            if ((time_start =
                 gsl_matrix_get(times_i, j, state_now)) != 0) {
                double dt = time_now - time_start;
                gsl_matrix_entry_add(hitting_times, j, state_now, dt);
                gsl_matrix_entry_add(times_n, j, state_now, 1.0);
                gsl_matrix_set(times_i, j, state_now, 0.0);
            }
            if ((jump_start =
                 gsl_matrix_get(jumps_i, j, state_now)) != 0) {
                unsigned int dj = jump_now - jump_start;
                gsl_matrix_entry_add(n_jumps,
                                     j, state_now, (double) dj);
                gsl_matrix_entry_add(jumps_n, j, state_now, 1.0);
                gsl_matrix_set(jumps_i, j, state_now, 0);
            }
        }
        c++;
        if (logl >= 1) {
            if (c % LOG_INTERVAL == 0 ||
                logl >= 3) {
                log_out << std::setiosflags(std::ios::fixed);
                log_out << "Analyze jump ";
                log_out << std::setprecision(0) << std::setw(8);
                log_out << c;
                log_out << "; time ";
                log_out << std::setprecision(1) << std::setw(7);
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
    // Now scale the invariant distribution and the matrizes.
    for (i = 0; i < ns; i++) {
        invariant_distribution[i] /=
            (p_times->back() - p_times->at(burn_in-1));
    }
    gsl_matrix_div_elements(hitting_times, times_n);
    gsl_matrix_div_elements(n_jumps, jumps_n);

    gsl_matrix_free(times_i);
}

void CTMCPath::print_hitting_times(std::ostream& out) {
    out << "Average hitting times." << std::endl;
    gsl_matrix_print(hitting_times, out);
}

void CTMCPath::print_n_jumps(std::ostream& out) {
    out << "Average number of jumps." << std::endl;
    gsl_matrix_print(n_jumps, out);
}

void CTMCPath::print_invariant_distribution(std::ostream& out) {
    out << "The invariant distribution is:" << std::endl;
    unsigned int i;
    out << std::setiosflags(std::ios::fixed);
    out << std::setprecision(5);
    for (i = 0; i < ns; i++) {
        out << std::setw(5) << invariant_distribution[i] << std::endl;
    }
}
