#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>
#include <boost/math/special_functions/binomial.hpp>

#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>

#include "ctmc.h"
#include "tools.h"
#include "getopt.h"

template double * new_zeroes<double>(unsigned int);

/// A K-dimensional state of the K-allelic Wright-Fisher process.
typedef std::vector<unsigned int> wfstate;
/// A flat (one dimensional) representation.
typedef unsigned int fstate;

void p_wfs(wfstate wfs, size_t K) {
    std::cout << "[ ";
    for (unsigned int i = 0; i < K; i++) {
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(0);
        std::cout << std::setw(2) << wfs.at(i);
        std::cout << " ";
    }
    std::cout << "]";
    std::cout << std::setprecision(10);
}

void init_states_three(wfstate fs_to_wfs[],
                       std::map<wfstate,fstate> &wfs_to_fs,
                       unsigned int N,
                       unsigned int S) {
    unsigned int c = 0;
    for (unsigned int i = 0; i <= N; i++) {
        for (unsigned int j = 0; j <= N; j++) {
            for (unsigned int k = 0; k <= N; k++) {
                if (i+j+k == N) {
                    wfstate wfs;
                    wfs.push_back(i);
                    wfs.push_back(j);
                    wfs.push_back(k);
                    fs_to_wfs[c] = wfs;
                    wfs_to_fs[wfs] = c;
                    c++;
                }
            }
        }
    }
    if (c != S) {
        out_error("Number of states does not match.");
    }
}

void init_states_four(wfstate fs_to_wfs[],
                      std::map<wfstate,fstate> &wfs_to_fs,
                      unsigned int N,
                      unsigned int S) {
    unsigned int c = 0;
    for (unsigned int i = 0; i <= N; i++) {
        for (unsigned int j = 0; j <= N; j++) {
            for (unsigned int k = 0; k <= N; k++) {
                for (unsigned int l = 0; l <= N; l++) {
                    if (i+j+k+l == N) {
                        wfstate wfs;
                        wfs.push_back(i);
                        wfs.push_back(j);
                        wfs.push_back(k);
                        wfs.push_back(l);
                        fs_to_wfs[c] = wfs;
                        wfs_to_fs[wfs] = c;
                        c++;
                    }
                }
            }
        }
    }
    if (c != S) {
        out_error("Number of states does not match.");
    }
}

// Set the diagonal elements so that row sum is sum.
void gsl_matrix_set_diag(gsl_matrix * m, double sum=0) {
    double row_sum;
    size_t s = m->size1;
    if (s != m->size2) out_error("Not a square matrix.");
    for (unsigned int i = 0; i < s; i++) {
        row_sum = 0;
        for (unsigned int j = 0; j < s; j++) {
            if (i != j) {
                row_sum += gsl_matrix_get(m, i, j);
            }
        }
        gsl_matrix_set(m, i, i, sum - row_sum);
    }    
}

gsl_matrix * params_to_mut_matrix_four(double m[], double pi[],
                                       double phi[], size_t N)
{
    gsl_matrix * u = gsl_matrix_alloc(4, 4);
    gsl_matrix_set(u, 0, 1, m[0]*pi[1] + phi[2]/pi[0]);
    gsl_matrix_set(u, 0, 2, m[1]*pi[2] - phi[1]/pi[0]);
    gsl_matrix_set(u, 0, 3, m[2]*pi[3] + phi[1]/pi[0] - phi[2]/pi[0]);
    gsl_matrix_set(u, 1, 0, m[0]*pi[0] - phi[2]/pi[1]);
    gsl_matrix_set(u, 1, 2, m[3]*pi[2] + phi[0]/pi[1]);
    gsl_matrix_set(u, 1, 3, m[4]*pi[3] - phi[0]/pi[1] + phi[2]/pi[1]);
    gsl_matrix_set(u, 2, 0, m[1]*pi[0] + phi[1]/pi[2]);
    gsl_matrix_set(u, 2, 1, m[3]*pi[1] - phi[0]/pi[2]);
    gsl_matrix_set(u, 2, 3, m[5]*pi[3] + phi[0]/pi[2] - phi[1]/pi[2]);
    gsl_matrix_set(u, 3, 0, m[2]*pi[0] - phi[1]/pi[3] + phi[2]/pi[3]);
    gsl_matrix_set(u, 3, 1, m[4]*pi[1] + phi[0]/pi[3] - phi[2]/pi[3]);
    gsl_matrix_set(u, 3, 2, m[5]*pi[2] - phi[0]/pi[3] + phi[1]/pi[3]);
    // gsl_matrix_print(u,std::cout,8);
    gsl_matrix_scale(u, (double) 1.0 / (double) N);
    gsl_matrix_set_diag(u, 1.0);
    // gsl_matrix_print(u,std::cout,8);
    return u;
}

gsl_matrix * print_rev_mut_matrix_four_params(double m[], double pi[],
                                              size_t N)
{
    gsl_matrix * u = gsl_matrix_alloc(4, 4);
    gsl_matrix_set(u, 0, 1, m[0]*pi[1]);
    gsl_matrix_set(u, 0, 2, m[1]*pi[2]);
    gsl_matrix_set(u, 0, 3, m[2]*pi[3]);
    gsl_matrix_set(u, 1, 0, m[0]*pi[0]);
    gsl_matrix_set(u, 1, 2, m[3]*pi[2]);
    gsl_matrix_set(u, 1, 3, m[4]*pi[3]);
    gsl_matrix_set(u, 2, 0, m[1]*pi[0]);
    gsl_matrix_set(u, 2, 1, m[3]*pi[1]);
    gsl_matrix_set(u, 2, 3, m[5]*pi[3]);
    gsl_matrix_set(u, 3, 0, m[2]*pi[0]);
    gsl_matrix_set(u, 3, 1, m[4]*pi[1]);
    gsl_matrix_set(u, 3, 2, m[5]*pi[2]);
    gsl_matrix_print(u,std::cout,8);
    gsl_matrix_scale(u, (double) 1.0 / (double) N);
    gsl_matrix_set_diag(u, 1.0);
    // gsl_matrix_print(u,std::cout,8);
    return u;
}

gsl_matrix * print_flux_param_matrix(double pi[],
                                     double phi[], size_t N)
{
    gsl_matrix * u = gsl_matrix_alloc(4, 4);
    gsl_matrix_set(u, 0, 1, + phi[2]/pi[0]/pi[1]);
    gsl_matrix_set(u, 0, 2, - phi[1]/pi[0]/pi[2]);
    gsl_matrix_set(u, 0, 3, (+ phi[1]/pi[0] - phi[2]/pi[0])/pi[3]);
    gsl_matrix_set(u, 1, 0, - phi[2]/pi[1]/pi[0]);
    gsl_matrix_set(u, 1, 2, + phi[0]/pi[1]/pi[2]);
    gsl_matrix_set(u, 1, 3, (- phi[0]/pi[1] + phi[2]/pi[1])/pi[3]);
    gsl_matrix_set(u, 2, 0, + phi[1]/pi[2]/pi[0]);
    gsl_matrix_set(u, 2, 1, - phi[0]/pi[2]/pi[1]);
    gsl_matrix_set(u, 2, 3, (+ phi[0]/pi[2] - phi[1]/pi[2])/pi[3]);
    gsl_matrix_set(u, 3, 0, (- phi[1]/pi[3] + phi[2]/pi[3])/pi[0]);
    gsl_matrix_set(u, 3, 1, (+ phi[0]/pi[3] - phi[2]/pi[3])/pi[1]);
    gsl_matrix_set(u, 3, 2, (- phi[0]/pi[3] + phi[1]/pi[3])/pi[2]);
    gsl_matrix_print(u,std::cout,8);
    gsl_matrix_scale(u, (double) 1.0 / (double) N);
    gsl_matrix_set_diag(u, 1.0);
    // gsl_matrix_print(u,std::cout,8);
    return u;
}

void wfs_to_x(wfstate wfs, size_t K, size_t N, gsl_vector * x) {
    for (unsigned int i = 0; i < K; i++) {
        double f = (double) wfs[i] / (double) N;
        gsl_vector_set(x,i,f);
    }
}

void get_x_new(gsl_vector * x_old, gsl_matrix * u,
                          gsl_vector * x_new) {
    gsl_blas_dgemv(CblasTrans, 1.0, u, x_old, 0.0, x_new);
}

void wfs_to_n(wfstate wfs, size_t K,
              unsigned int n[]) {
    for (unsigned int i = 0; i < K; i++) {
        n[i] = wfs.at(i);
    }
}

double transition_prob(fstate i, fstate j,
                       wfstate fs_to_wfs[],
                       gsl_matrix * u,
                       size_t K, size_t N) {
    gsl_vector * x_old = gsl_vector_alloc(K);
    gsl_vector * x_new = gsl_vector_alloc(K);
    wfstate wfsa = fs_to_wfs[i];
    wfs_to_x(wfsa, K, N, x_old);
    get_x_new(x_old, u, x_new);
    wfstate wfsb = fs_to_wfs[j];
    unsigned int n[K];
    wfs_to_n(wfsb, K, n);
    double p = gsl_ran_multinomial_pdf(K, x_new->data, n);
    
    // std::cout << "Transition probability from ";
    // p_wfs(wfsa, K);
    // std::cout << " to ";
    // p_wfs(wfsb, K);
    // std::cout << std::setprecision(10);
    // std::cout << ": " << p << std::endl;
    
    gsl_vector_free(x_old);
    gsl_vector_free(x_new);
    return p;
}

gsl_matrix * general_wright_fisher_mut_matrix
(wfstate fs_to_wfs[],
 gsl_matrix * u,
 size_t K,
 size_t N,
 unsigned int S)
{
    gsl_matrix * q = gsl_matrix_alloc(S, S);
    gsl_matrix_set_zero(q);

    // Set the rates of frequency shifts.
    for (fstate i = 0; i < S; i++) {
        for (fstate j = 0; j < S; j++) {
            if (i != j) {
                gsl_matrix_set
                    (q, i, j,
                     transition_prob(i, j, fs_to_wfs, u, K, N));
            }
        }
    }
    gsl_matrix_set_diag(q);
    return q;
}

double trans_rate(gsl_matrix * q, wfstate wfsa, wfstate wfsb,
                  std::map<wfstate,fstate> wfs_to_fs) {
    fstate i = wfs_to_fs[wfsa];
    fstate j = wfs_to_fs[wfsb];
    return gsl_matrix_get(q, i, j);
}

int main()
{
    //////////////////////////////
    /// The number of alleles.  Only works for K=3 at the moment.
    size_t K = 4;
    /// The population size.
    size_t N = 30;
    /// Simulation time.
    double tm = 5e9;
    /// Parameters (mu_ab, pi_a, 0.5*phi_ab); cf. Burden and Tang (2016),
    /// 4 allele case.
    double m[] = {1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2};
    double pi[] = {0.1, 0.2, 0.3, 0.4};
    double phih[] = {0.2e-2, 0.05e-2, -0.015e-2};
    /// The log level.
    int log_l = 1;
    /// Set to true for random seed.
    bool seed = false;

    //////////////////////////////
    // The total number of states.
    double Sd = boost::math::binomial_coefficient<double>(N+K-1, N);
    unsigned int S = (unsigned int) Sd;
    // It is natural, to have a K-dimensional representation of a
    // state of the K-allelic Wright-Fisher process.  However, the
    // mutation matrix is still two-dimensional (flat).  These
    // variables convert a flat one-dimensional state to a
    // multidimensional WF-state and vice versa.
    wfstate fs_to_wfs[S];
    std::map<wfstate,fstate> wfs_to_fs;
    if (K == 3) init_states_three(fs_to_wfs, wfs_to_fs, N, S);
    else if (K == 4) init_states_four(fs_to_wfs, wfs_to_fs, N, S);
    else out_error("Only three and four alleles supported.");
    // TODO: Write this function for K=3.
    gsl_matrix * u = NULL;
    if (K == 3) out_error("NOT IMPLEMENTED.");
    else if (K == 4)
        u = params_to_mut_matrix_four(m, pi, phih, N);
    else out_error("Only three and four alleles supported.");

    print_flux_param_matrix(pi, phih, N);
    exit(1);
    
    std::cout << "Compute transition rate matrix." << std::endl;
    gsl_matrix * q =
        general_wright_fisher_mut_matrix(fs_to_wfs, u, K, N, S);
    CTMC chain(q, S);
    chain.set_log_level(log_l);
    // chain.print_info(std::cout);
    if (seed) {
        unsigned long int s;
        std::cout << "Bytes set:";
        std::cout <<
            syscall(SYS_getrandom, &s, sizeof(unsigned long int), 0);
        std::cout << std::endl;
        std::cout << "Seed is: " << s << "." << std::endl;
        chain.rg->set_seed(s);
    }
    chain.burn_it_in();

    std::cout << "Run chain." << std::endl;
    chain.run(tm);

    std::cout << "Print output." << std::endl;
    // chain.print_direct_hitting_times(std::cout);
    // chain.print_hitting_times(std::cout);
    // chain.print_direct_number_jumps(std::cout);
    // chain.print_invariant_distribution(std::cout);

    std::cout << "A1-A2 edge." << std::endl;
    for (unsigned int i = 0; i <= N; i++) {
        wfstate wfs;
        wfs.resize(K);
        wfs[0] = i;
        wfs[1] = N-i;
        wfs[2] = 0;
        wfs[3] = 0;
        p_wfs (wfs, K);
        std::cout << " ";
        std::cout << std::log10
            (N * chain.get_entry_invariant_distribution(wfs_to_fs[wfs]));
        std::cout << std::endl;
    }
    std::cout << "A1-A3 edge." << std::endl;
    for (unsigned int i = 0; i <= N; i++) {
        wfstate wfs;
        wfs.resize(K);
        wfs[0] = i;
        wfs[1] = 0;
        wfs[2] = N-i;
        wfs[3] = 0;
        p_wfs (wfs, K);
        std::cout << " ";
        std::cout << std::log10
            (N * chain.get_entry_invariant_distribution(wfs_to_fs[wfs]));
        std::cout << std::endl;
    }
    std::cout << "A1-A4 edge." << std::endl;
    for (unsigned int i = 0; i <= N; i++) {
        wfstate wfs;
        wfs.resize(K);
        wfs[0] = i;
        wfs[1] = 0;
        wfs[2] = 0;
        wfs[3] = N-i;
        p_wfs (wfs, K);
        std::cout << " ";
        std::cout << std::log10
            (N * chain.get_entry_invariant_distribution(wfs_to_fs[wfs]));
        std::cout << std::endl;
    }
    std::cout << "A2-A3 edge." << std::endl;
    for (unsigned int i = 0; i <= N; i++) {
        wfstate wfs;
        wfs.resize(K);
        wfs[0] = 0;
        wfs[1] = i;
        wfs[2] = N-i;
        wfs[3] = 0;
        p_wfs (wfs, K);
        std::cout << " ";
        std::cout << std::log10
            (N * chain.get_entry_invariant_distribution(wfs_to_fs[wfs]));
        std::cout << std::endl;
    }
    std::cout << "A2-A4 edge." << std::endl;
    for (unsigned int i = 0; i <= N; i++) {
        wfstate wfs;
        wfs.resize(K);
        wfs[0] = 0;
        wfs[1] = i;
        wfs[2] = 0;
        wfs[3] = N-i;
        p_wfs (wfs, K);
        std::cout << " ";
        std::cout << std::log10
            (N * chain.get_entry_invariant_distribution(wfs_to_fs[wfs]));
        std::cout << std::endl;
    }
    std::cout << "A3-A4 edge." << std::endl;
    for (unsigned int i = 0; i <= N; i++) {
        wfstate wfs;
        wfs.resize(K);
        wfs[0] = 0;
        wfs[1] = 0;
        wfs[2] = i;
        wfs[3] = N-i;
        p_wfs (wfs, K);
        std::cout << " ";
        std::cout << std::log10
            (N * chain.get_entry_invariant_distribution(wfs_to_fs[wfs]));
        std::cout << std::endl;
    }

    gsl_matrix_free(u);
    gsl_matrix_free(q);
    return 0;
}
