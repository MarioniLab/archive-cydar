#include "cydar.h"
#include "utils.h"

matrix_info::matrix_info (int nr, int nc, const double* ptr) : nrow(nr), ncol(nc), dptr(ptr) {}

matrix_info check_matrix(SEXP matrix) {
    if (!isReal(matrix)) {
        throw std::runtime_error("matrix must be double-precision");
    }

    SEXP dims=getAttrib(matrix, R_DimSymbol);
    if (!isInteger(dims) || LENGTH(dims)!=2) { 
        throw std::runtime_error("dimensions of the matrix should be an integer vector of length 2");
    }
    int nrow=INTEGER(dims)[0], ncol=INTEGER(dims)[1];
    if (LENGTH(matrix)!=nrow*ncol) {
        throw std::runtime_error("recorded dimensions of the matrix are not consistent with its length"); 
    }

    matrix_info output(nrow, ncol, REAL(matrix));
    return output;
}


