#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/** 
 * Print a warning to std::cerr.
 * 
 * @param warn warning message.
 */
void out_warning(const char * warn);

/** 
 * Print a warning to std::cerr.
 * 
 * @param warn warning message.
 */
void out_warning(std::string warn);

/** 
 * Print an error to std::cerr and exit.
 * 
 * @param error error message.
 */
void out_error(const char * error);

/** 
 * Print an error to std::cerr and exit.
 * 
 * @param error error message.
 */
void out_error(std::string error);

/** 
 * Initialize an std::vector with zeroes.
 * 
 * @param length length of the vector. 
 * 
 * @return the vector.
 */
template<typename T>
T * new_zeroes(unsigned int length) {
    T * v = new T[length];
    for (unsigned int i = 0; i < length; i++) {
        v[i] = 0;
    }
    return v;
}

/** 
 * Set a row of a gsl_matrix to a value.
 * 
 * @param m pointer to matrix.
 * @param i row number.
 * @param x new value.
 */
void gsl_matrix_set_row_to_value (gsl_matrix * m,
                                  const size_t i,
                                  double x);


/** 
 * Add value to element in matrix.
 * 
 * @param m pointer to matrix.
 * @param i row number.
 * @param j column number.
 * @param x value to add.
 */
void gsl_matrix_entry_add(gsl_matrix * m,
                          const size_t i,
                          const size_t j,
                          double x);

/** 
 * Print a gsl_matrix.
 * 
 * @param m the matrix to print.
 */
void gsl_matrix_print(gsl_matrix * m, std::ostream& out,
                      size_t precision=3);

#endif
