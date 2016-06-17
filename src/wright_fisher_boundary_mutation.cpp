#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <string>
#include "ctmc.h"
#include "tools.h"
#include "getopt.h"
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>

template double * new_zeroes<double>(unsigned int);



gsl_matrix * wright_fisher_mut_matrix (size_t n, double mu) {
    gsl_matrix * m = gsl_matrix_alloc(n+1, n+1);
    gsl_matrix_set_zero(m);

    // Set the rates of frequency shifts.
    unsigned int i;
    unsigned int j;
    //////////////////////////////
    // Boundary mutation only.
    // Starting at state i; the frequency of the first allele is i/n.
    // Going to state j; picking j times the first allele.
    // Starting at 0 and n is not necessary, because then only mutations can happen.
    for (i = 1; i < n; i++) {
        for (j = 0; j <= n; j++) {
            double p = (double) i / (double) n;
            gsl_matrix_set(m, i, j, gsl_ran_binomial_pdf (j, p, n));
        }
    }
    // Set the mutations from the boundaries.
    for (j = 0; j <= n; j++) {
        double smu = mu / (double) n;
        gsl_matrix_set(m, 0, j, gsl_ran_binomial_pdf (j, smu, n));
        gsl_matrix_set(m, n, j, gsl_ran_binomial_pdf (j, (1-smu), n));
    }

    // //////////////////////////////
    // // General mutations.
    // // Starting at state i; the frequency of the first allele is i/n.
    // // Going to state j; picking j times the first allele.
    // for (i = 0; i <= n; i++) {
    //     double f = (double) i / (double) n;
    //     double smu = mu / (double) n;
    //     double p = ( f*(1-smu) + (1-f)*smu);
    //     for (j = 0; j <= n; j++) {
    //         gsl_matrix_set(m, i, j, gsl_ran_binomial_pdf (j, p, n));
    //     }
    // }

    //////////////////////////////
    // Set the diagonal elements.
    double row_sum;
    for (i = 0; i <= n; i++) {
        row_sum = 0;
        for (j = 0; j <= n; j++) {
            if (i != j) {
                row_sum += gsl_matrix_get(m,i,j);
            }
        }
        gsl_matrix_set(m,i,i,-row_sum);
    }

    return m;
}

int main(int argc, char *argv[])
{
    // Option parsing.
    unsigned int ne = 10;
    int tm = 1e6;
    double mu = 1e-2;
    int log_l = 1;
    char * out_fn = NULL;
    bool seed = false;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "n:t:m:f:l:s")) != -1)
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
            case 's':
                // Set seed randomly.
                seed = true;
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

    log_out << "Setup chain." << std::endl;
    gsl_matrix * m =
        wright_fisher_mut_matrix(ne, mu);
    CTMC chain(m, ne+1);
    if (seed) {
        unsigned long int s;
        log_out << "Bytes set:";
        log_out <<
            syscall(SYS_getrandom, &s, sizeof(unsigned long int), 0);
        log_out << std::endl;
        log_out << "Seed is: " << s << "." << std::endl;
        chain.rg->set_seed(s);
    }
    chain.burn_it_in();
    chain.set_log_level(log_l);
    // chain.print_info(std::cout);

    std::cout << "Run chain." << std::endl;
    chain.run(tm);

    std::cout << "Print output." << std::endl;
    // chain.print_direct_hitting_times(std::cout);
    chain.print_hitting_times(log_out);
    // chain.print_direct_number_jumps(std::cout);
    chain.print_invariant_distribution(log_out);

    gsl_matrix_free(m);
    log_out.close();
    return 0;
}
