#ifndef UTILS_H
#define UTILS_H

// Check matrix inputs.

struct matrix_info {
    matrix_info(int, int, bool);
    const size_t nrow, ncol;
    const bool is_integer;
    const int* iptr;
    const double* dptr;
};

matrix_info check_matrix(SEXP matrix);

#endif
