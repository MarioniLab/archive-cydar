#include "cydar.h"
#include "objects.h"

/*****************************************
 * reassign the cells that were counted back to the hypersphere centre.
 *****************************************/

SEXP find_counted(SEXP hyperx, SEXP centers, SEXP info, SEXP cell_exprs, SEXP distance) try {
    finder fx(hyperx, R_NilValue, centers, info);
    matrix_info cx=check_matrix(cell_exprs); 
    const size_t& ncells=cx.ncol;
    const size_t& nmarkers=cx.nrow;
    if (nmarkers!=fx.searcher->get_nmarkers()) { 
        throw std::runtime_error("number of markers is not consistent between matrices");
    }
    if (!isReal(distance)|| LENGTH(distance)!=1) { throw std::runtime_error("distance must be a double-precision scalar"); }
    const double threshold=asReal(distance);

    SEXP overall_output=PROTECT(allocVector(LGLSXP, ncells));
    try {
        int* optr=LOGICAL(overall_output);
        const double* curexprs=cx.dptr;
        std::deque<double>& distances=fx.searcher->distances;

        for (size_t c=0; c<ncells; ++c, curexprs+=nmarkers) {
            fx.searcher->find_nearest_neighbors(curexprs, 1, true); 
            if (distances.back() <= threshold) { 
                optr[c]=1;
            } else {
                optr[c]=0;
            }
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);
    return overall_output;
} catch (std::exception& e) {
    return mkString(e.what());
}

