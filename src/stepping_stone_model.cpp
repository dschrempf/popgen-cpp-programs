/**
 * @file   stepping_stone_model.cpp
 * @author Dominik Schrempf <dominik.schrempf@gmail.com>
 * @date   Fri Mar 13 16:56:17 2015
 *
 * @brief  Simulate Stepping Stone Model with a Markov chain
 *
 * Stepping Stone Model: In this model we have an n-by-n array of
 *  squares, and each square is initially any one of k different
 *  colors. For each step, a square is chosen at random. This square
 *  then chooses one of its eight neighbors at random and assumes the
 *  color of that neighbor.
 *
 *  To avoid boundary problems, we assume that if a square S is on the
 *  left-hand boundary, say, but not at a corner, it is adjacent to
 *  the square T on the right-hand boundary in the same row as S, and
 *  S is also adjacent to the squares just above and below T . A
 *  similar assumption is made about squares on the upper and lower
 *  boundaries. (These adjacencies are much easier to understand if
 *  one imagines making the array into a cylinder by gluing the top
 *  and bottom edge together, and then making the cylinder into a
 *  doughnut by gluing the two circular boundaries together.) With
 *  these adjacencies, each square in the array is adjacent to exactly
 *  eight other squares.
 *
 *  A state in this Markov chain is a description of the color of each
 *  square. For this 2 Markov chain the number of states is k n ,
 *  which for even a small array of squares is enormous. This is an
 *  example of a Markov chain that is easy to simulate but difficult
 *  to analyze in terms of its transition matrix.
 *
 */

#include <iostream>
#include <stdexcept>
#include <gsl/gsl_rng.h>

// User interface.
int l = 30;                    /**< Array length; there are n*n squares. */

// Markov chain settings.
int n_steps = 100000;

void
print_field (char * field, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n ; j++) {
            if (field[i*n+j] == 1) std::cout << "+";
            else if (field[i*n+j] == 0) std::cout << " ";
            else throw std::runtime_error("Unknown field value.\n");
        }
        std::cout << std::endl;
    }
}

void jump(char * field, int n, const gsl_rng * r)
{
    int nn = n*n;
    volatile int pick = gsl_rng_uniform_int (r, nn);
    volatile int direction = gsl_rng_uniform_int (r, 8);

    volatile int i = (int) pick / n;
    volatile int j = pick % n;

    switch(direction) {
    case 0 :
        i -= 1;
        j -= 1;
        break;
    case 1 :
        j -= 1;
        break;
    case 2 :
        i += 1;
        j -= 1;
        break;
    case 3 :
        i += 1;
        break;
    case 4 :
        i += 1;
        j += 1;
        break;
    case 5 :
        j += 1;
        break;
    case 6 :
        i -= 1;
        j += 1;
        break;
    case 7 :
        i -= 1;
        break;
    default :
        throw std::runtime_error("Direction not known.");
    }
    if (i < 0) i += n;
    else if (i >= n) i -= n;

    if (j < 0) j += n;
    else if (j >=n ) j -= n;

    field[pick] = field[i*n+j];
}

int
main (void)
{
    const gsl_rng_type * t;
    gsl_rng * r;

    gsl_rng_env_setup();
    t = gsl_rng_default;
    r = gsl_rng_alloc (t);
    gsl_rng_set (r, 1);

    int ll = l*l;
    char field[ll];
    for (int i = 0; i < ll; i++) {
        field[i] = i % 2;
    }

    for (int i = 0; i < n_steps; i++) {
        jump (field, l, r);
    }


    print_field (field, l);

    gsl_rng_free (r);

    return 0;
}
