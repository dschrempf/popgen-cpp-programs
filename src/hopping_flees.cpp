#include <iostream>
#include <iomanip>
#include <cmath>
#include "ran_generator.h"

/**
 * @file   hopping_flees.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Sep 23 11:47:40 2015
 *
 * @brief  Markov Chains - Norris; exercise 2.8.1; Flees hop around a triangle.
 *
 * Two fleas are bound together to take part in a nine-legged race on
 * the vertices A, B, C of a triangle.  Feal 1 hops at random times in
 * the clockwise direction; each hop takes the pair from one vertex to
 * the next and the times between successive hops of Flea 1 are
 * independent random variables, each with exponential distribution,
 * mean \f$ 1/\lambda \f$.  Flea 2 behaves similarly, but hops in the
 * anticlockwise direction, the times between his hops having mean \f$
 * 1/\mu \f$.
 *
 *
\verbatim
      C
   
  A       B
\endverbatim
 * 
 * We can also calculate the analytical solution for this problem
 * which is of the form
 * \f{align}{
 *  \frac{1}{3} + \frac{2}{3} \text{exp}\left(-\frac{3(\lambda + \mu)t}{2}\right)
 *  \text{cos}\left( \frac{\sqrt{3}(\lambda - \mu)t}{2}\right).
 * \f}
 *
 * We can see, that the empirical simulation fits the analytical results very well:
 *
 * \image html hopping_flees.png
 * 
 */

size_t n_states     = 3;
unsigned int s_init = 0;

double mu = 1.0;
double la = 2.0;

double q = mu+la;
double q_matrix[] = { -q, mu, la,
                      la, -q, mu,
                      mu, la, -q};

double t_max              = 1;
unsigned int n_time_steps = 1000;
unsigned int n_reps       = 10000;

class ContChain
{
public:
    ContChain (int s_init, double t_init, int d_init, double * qm_init);
    ~ContChain ();
    int s;
    double t;
    int d;
    double * qm;
};

ContChain::ContChain (int s_init, double t_init, int d_init, double * qm_init)
{
    s = s_init;
    t = t_init;
    d = d_init;
    qm = qm_init;
}

ContChain::~ContChain ()
{

}

bool jump(ContChain * ch, double t, RanGen * rg)
{
    int s = ch->s;
    int d = ch->d;
    double q = -ch->qm[s*d + s];
    // std::cout << "q=" << q << std::endl;
    double dt = rg->pick_exponential(1.0/q);
    // std::cout << "dt=" << dt << std::endl;
    if (ch->t + dt > t) {
        ch->t = t;
        // std::cout << "Stop at " << t << std::endl;
        return false;
    }
    else {
        // FIXME: Only works for this specific problem at the moment.
        gsl_vector * v = gsl_vector_alloc(2);
        gsl_vector_set(v,0,la);
        gsl_vector_set(v,1,mu);
        // v[] = {la, mu};
        int ds = rg->vector_weighted_pick(v, 2);
        gsl_vector_free(v);
        if (ds == 0) ds = -1;
        int s_prime = s + ds;
        if (s_prime < 0) s_prime = 2;
        if (s_prime > 2) s_prime = 0;
        ch->s  = s_prime;
        ch->t += dt;
        return true;
    }
}

void evolve_chain(ContChain * ch, double t, RanGen * rg)
{
    while (jump(ch,t,rg));
}

double hopping_flees_prob (double la, double mu, double t) {
    double e = exp (-3.0*(la+mu)*t/2.0);
    double c = cos (sqrt(3.0)*(la-mu)*t/2.0);
    return 1.0/3.0 + 2.0/3.0*e*c;
}

int main()
{
    RanGen rg;

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    std::cout << std::setw(7) << "t";
    std::cout << std::setw(7) << "p_emp";
    std::cout << std::setw(7) << "p_exa" << std::endl;

    unsigned int hits;
    double p_emp;
    double p_exa;

    for (unsigned int i = 0; i < n_time_steps; i++) {
        hits = 0;
        double t = (double) i * t_max / n_time_steps;
        for (unsigned int r = 0; r < n_reps; r++) {
            ContChain * ch = new ContChain(s_init, 0.0, n_states, q_matrix);
            evolve_chain (ch, t, &rg);
            if (ch->s == 0) hits++;
            delete ch;
        }
        p_emp = (double) hits/n_reps;
        p_exa = hopping_flees_prob(la, mu, t);

        std::cout << std::setw(7) << t;
        std::cout << std::setw(7) << p_emp;
        std::cout << std::setw(7) << p_exa << std::endl;
    }

    return 0;
}
