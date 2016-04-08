#include <gsl/gsl_matrix.h>
#include <iostream>
#include <string>
#include "ctmc.h"
#include "tools.h"
#include "getopt.h"

template double * new_zeroes<double>(unsigned int);

gsl_matrix * moran_boundary_mut_matrix (size_t n, double mu) {
    gsl_matrix * m = gsl_matrix_alloc(n+1, n+1);
    gsl_matrix_set_zero(m);

    // Set the rates of frequency shifts.
    size_t i;
    for (i = 1; i < n; i++) {
        gsl_matrix_set(m, i, i-1, (double) i*(n-i)/n);
        gsl_matrix_set(m, i, i, (double) -2*i*(n-i)/n);
        gsl_matrix_set(m, i, i+1, (double) i*(n-i)/n);
    }

    // Set the mutations from the boundaries.
    gsl_matrix_set(m, 0, 0, -mu);
    gsl_matrix_set(m, 0, 1, mu);
    gsl_matrix_set(m, n, n, -mu);
    gsl_matrix_set(m, n, n-1, mu);

    return m;
}

int main(int argc, char *argv[])
{
    // Option parsing.
    unsigned int ne = 10;
    int tm = 1e6;
    double mu = 1e-1;
    int log_l = 1;
    char * out_fn = NULL;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "n:t:m:f:l:")) != -1)
        switch (c)
            {
            case 'n':
                // Population size.
                ne = atoi(optarg);
                break;
            case 't':
                // Run time of CTMC.
                tm = atoi(optarg);
                break;
            case 'm':
                // Mutation rate.
                mu = atof(optarg);
                break;
            case 'f':
                // Print output to file.
                out_fn = optarg;
                break;
            case 'l':
                // Set log level.
                log_l = atoi(optarg);
                break;
            case '?':
                if (optopt == 'n' || optopt == 't' || optopt == 'm') {
                    std::cerr << "Option -" << optopt;
                    std::cerr << " requires an argument." << std::endl;
                }
                else {
                    std::cerr << "Unknown option `-" << optopt;
                    std::cerr << "'.\n" << std::endl;
                }
                return 1;
            default:
                abort ();
            }

    // Open log if given.
    if (out_fn == NULL)
        out_error("Please give a log file with -f LOG_FILE.");

    std::ofstream log_out;
    log_out.open(out_fn, std::ios::out);

    // Print command line.
    log_out << "Command line: ";
    for (int count=0; count < argc; ++count)
        log_out << argv[count] << " ";
    log_out << std::endl;

    gsl_matrix * m =
        moran_boundary_mut_matrix(ne, mu);
    CTMC chain(m, ne+1);
    // CTMC chain(m, ne+1, true);
    chain.set_log_level(log_l);
    // chain.print_info(std::cout);

    chain.run(tm);
    // chain.print_direct_hitting_times(std::cout);
    chain.print_hitting_times(log_out);
    // chain.print_direct_number_jumps(std::cout);
    chain.print_invariant_distribution(log_out);

    gsl_matrix_free(m);
    log_out.close();
    return 0;
}
