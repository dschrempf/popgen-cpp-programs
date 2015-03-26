/**
 * @file   hitchhiking.c
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Wed Feb 11 14:31:40 2015
 *
 * @brief  Simulate hitchhiking along a positively selected locus.
 *
 * Gillespie - Population Genetics
 * p. 110
 * Problem 4.3
 *
 * Write a computer simulation to simulate the hitchhiking along a
 * positively selected locus and look at the heterozygousity near this
 * locus.
 *
 * The parameters are:
 * - \f$p_a\f$ frequency of \f$A1\f$
 * - \f$q_a\f$ frequency of \f$A2\f$
 * - \f$p_b\f$ frequency of \f$B1\f$
 * - \f$q_b\f$ frequency of \f$B2\f$
 *
 * - \f$x_1\f$ frequency of \f$A_1B_1\f$
 * - \f$x_2\f$ frequency of \f$A_1B_2\f$
 * - \f$x_3\f$ frequency of \f$A_2B_1\f$
 * - \f$x_4\f$ frequency of \f$A_2B_2\f$

 * The linkage disequlibrium is:
 * - \f$D = x_1 - p_A * p_B\f$
 *
 */



#include <stdio.h>

#define h (0.5)		/* heterozygous coefficient */
#define s (0.1)		/* selection coefficient */

/** 
 * Marginal fitness 1.
 * 
 * @param p 
 * 
 * @return 
 */
double w_bar_1 (double p)
{
    double q = 1 - p;
    return 1.0 - q*h*s;
}

/** 
 * Marginal fitness 3.
 * 
 * @param p 
 * 
 * @return 
 */
double w_bar_3 (double p)
{
    double q = 1 -p;
    return 1.0 - p*h*s - q*s;
}

/** 
 * Mean fitness.
 * 
 * @param p 
 * 
 * @return 
 */
double w_bar (double p)
{
  double q = 1-p;
  return 1 - 2*p*q*h*s - q*q*s;
}

/** 
 * Change in frequency of A1B1.
 * 
 * @param x1 
 * @param x2 
 * @param x3 
 * @param x4 
 * @param r 
 * 
 * @return 
 */
double delta_x1 (double x1, double x2, double x3, double x4, double r) {
  double p = x1 + x2;
  double q = 1-p;
  double num;

  num = x1 * (w_bar_1 (p) - w_bar (p)) - r*(1-h*s) * (x1*x4 - x2*x3);

  return num / w_bar (p);
}

/** 
 * Change in frequency of A1B2.
 * 
 * @param x1 
 * @param x2 
 * @param x3 
 * @param x4 
 * @param r 
 * 
 * @return 
 */
double delta_x2 (double x1, double x2, double x3, double x4, double r) {
  double p = x1 + x2;
  double q = 1-p;
  double num;

  num = x2 * (w_bar_1 (p) - w_bar (p)) + r*(1-h*s) * (x1*x4 - x2*x3);

  return num / w_bar (p);
}

/** 
 * Change in frequency of A2B1.
 * 
 * @param x1 
 * @param x2 
 * @param x3 
 * @param x4 
 * @param r 
 * 
 * @return 
 */
double delta_x3 (double x1, double x2, double x3, double x4, double r) {
  double p = x1 + x2;
  double q = 1-p;
  double num;

  num = x3 * (w_bar_3 (p) - w_bar (p)) + r*(1-h*s) * (x1*x4 - x2*x3);

  return num / w_bar (p);
}

/** 
 * Change in frequency of A2B2.
 * 
 * @param x1 
 * @param x2 
 * @param x3 
 * @param x4 
 * @param r 
 * 
 * @return 
 */
double delta_x4 (double x1, double x2, double x3, double x4, double r) {
  double p = x1 + x2;
  double q = 1-p;
  double num;

  num = x4 * (w_bar_3 (p) - w_bar (p)) - r*(1-h*s) * (x1*x4 - x2*x3);

  return num / w_bar (p);
}

/** 
 * Iterate the equations and return the final heterozygosity.
 * 
 * @param x1 
 * @param x2 
 * @param x3 
 * @param x4 
 * @param r 
 * @param eps 
 * 
 * @return 
 */
double get_het_f (double x1, double x2,
		  double x3, double x4,
		  double r, double eps) {
  double pa;
  double qa;
  double pb;
  double qb;
  double het_f;

  double dx1;
  double dx2;
  double dx3;
  double dx4;

  while (x1+x2 < 1.0 - eps) {
    pa = x1 + x2;
    qa = 1.0 - pa;

    dx1 = delta_x1 (x1, x2, x3, x4, r);
    dx2 = delta_x2 (x1, x2, x3, x4, r);
    dx3 = delta_x3 (x1, x2, x3, x4, r);
    dx4 = delta_x4 (x1, x2, x3, x4, r);

    x1 = x1 + dx1;
    x2 = x2 + dx2;
    x3 = x3 + dx3;
    x4 = x4 + dx4;
  }
  pb = x1 + x3;
  qb = 1.0 - pb;
  het_f = 2.0 * pb * qb;

  return het_f;
}

int main(void) {
  int i;
  FILE *fp;

  /* Diploid population size. */
  int N = 500;
  /* Recombination rate. */
  double r;
  /* Threshold used to start and stop iteration. */
  double eps = 1.0 / (2.0 * (double) N);

  /* Genotype frequencies. */
  double x1 = eps;
  double x2 = 0;
  double x3 = 0.5 - x1;
  double x4 = 0.5;

  double het_i;
  double het_f;

  het_i = 2.0 * (x1 + x3) * (1 - x1 - x3);

  fp = fopen("hitchhikingN500.dat", "w");

  for (i = 0; i < 100; i++) {
    r = (double) i / 2000;
    het_f = get_het_f (x1, x2, x3, x4, r, eps);
    fprintf(fp, "%f\t%f\n", r/s, het_f / het_i);
  }

  return 0;
}
