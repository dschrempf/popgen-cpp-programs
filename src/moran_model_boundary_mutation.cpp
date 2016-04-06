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
    int ne = 10;
    int tm = 1e6;
    double mu = 1e-1;
    char * ht_fn = NULL;
    char * nj_fn = NULL;
    char * id_fn = NULL;
    int logl = 1;
    int c;
    
    opterr = 0;

    while ((c = getopt (argc, argv, "n:t:m:h:j:i:l:")) != -1)
        switch (c)
            {
            case 'n':
                ne = atoi(optarg);
                break;
            case 't':
                tm = atoi(optarg);
                break;
            case 'm':
                mu = atof(optarg);
                break;
            case 'h':
                ht_fn = optarg;
                break;
            case 'j':
                nj_fn = optarg;
                break;
            case 'i':
                id_fn = optarg;
                break;
            case 'l':
                logl = atoi(optarg);
                break;
            case '?':
                if (optopt == 'n' || optopt == 't' || optopt == 'm')
                    std::cerr << "Option -" << optopt << " requires an argument." << std::endl;
                else
                    std::cerr << "Unknown option `-" << optopt << "'.\n" << std::endl;
                return 1;
            default:
                abort ();
            }
    
    // Print command line.
    for (int count=0; count < argc; ++count)
        std::cout << argv[count] << " ";
    std::cout << std::endl;
    
    gsl_matrix * m =
        moran_boundary_mut_matrix(ne, mu);
    CTMC chain(m, ne+1, std::cout);
    chain.set_log_level(logl);
    // chain.activate_logs();
    // chain.print_info(std::cout);

    chain.evolve(tm);
    // chain.pretty_print_log();
    CTMCPath * log = chain.ctmc_path;
    log->analyze();

    if (ht_fn != NULL) {
        std::ofstream ht_of;
        ht_of.open(ht_fn, std::ios::out);
        log->print_hitting_times(ht_of);
        ht_of.close();
    }
    else log->print_hitting_times(std::cout);
    if (nj_fn != NULL) {
        std::ofstream nj_of;
        nj_of.open(nj_fn, std::ios::out);
        log->print_n_jumps(nj_of);
        nj_of.close();
    }
    else log->print_n_jumps(std::cout);
    if (id_fn != NULL) {
        std::ofstream id_of;
        id_of.open(id_fn, std::ios::out);
        log->print_invariant_distribution(id_of);
        id_of.close();
    }
    else log->print_invariant_distribution(std::cout);
                       
    gsl_matrix_free(m);
    return 0;
}
