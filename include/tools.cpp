#include "tools.h"

void out_warning(std::string warn) {
    out_warning(warn.c_str());
}

void out_warning(const char * warn) {
    std::cerr << "WARNING: " << warn << std::endl;
}

void gsl_matrix_set_row_to_value(gsl_matrix * m,
                                 const size_t i,
                                 double x) {
    size_t row_length = m->size2;
    gsl_vector * v = gsl_vector_alloc(row_length);
    gsl_vector_set_all(v, x);
    gsl_matrix_set_row(m, i, v);
    gsl_vector_free(v);
}

void gsl_matrix_entry_add(gsl_matrix * m,
                          const size_t i,
                          const size_t j,
                          double x) {
    double old = gsl_matrix_get(m, i, j);
    gsl_matrix_set(m, i, j, old+x);
}

void gsl_matrix_print(gsl_matrix * m, std::ostream& out,
                      size_t p) {
    unsigned int s1 = m->size1;
    unsigned int s2 = m->size2;
    unsigned int i,j;
    out << std::setiosflags(std::ios::fixed);
    out << std::setprecision(p);
    for (i = 0; i < s1; i++) {
        for (j = 0; j < s2; j++) {
            out << std::setw(p+6) << gsl_matrix_get(m, i, j);
        }
        out << std::endl;
    }
}

void f_ostream_to_file(void (*p_f)(std::ostream&),
                       char const * fn) {
    std::ofstream outfile;
    outfile.open(fn, std::ios::out);
    p_f(outfile);
    outfile.close();
}
