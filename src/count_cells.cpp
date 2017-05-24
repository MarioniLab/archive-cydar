#include "objects.h"
#include "packer.h"

SEXP count_cells(SEXP exprs, SEXP distance, SEXP centers, SEXP cluster_info, SEXP curcells) try {
    auto searcher=generate_holder(exprs, centers, cluster_info);
    std::deque<size_t>& collected=searcher->neighbors;

    // Checking distances and chosen cells.
    if (!isReal(distance)|| LENGTH(distance)!=1) { throw std::runtime_error("distance must be a double-precision scalar"); }
    const double threshold=asReal(distance);
    if (!isInteger(curcells)) {
        throw std::runtime_error("chosen indices should be an integer vector");
    }
    const int N=LENGTH(curcells);
    const int ncells=searcher->get_ncells();
    const size_t& nmarkers=searcher->get_nmarkers();
    if (nmarkers==0) {
        throw std::runtime_error("number of markers should be positive");
    }
    const int* cptr=INTEGER(curcells);
   
    // Setting up output vectors. 
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, N));
        SEXP cellids=VECTOR_ELT(output, 0);
        std::deque<int> sorted_ids;

        SET_VECTOR_ELT(output, 1, allocVector(INTSXP, N));
        int* numptr=INTEGER(VECTOR_ELT(output, 1));

        // Running through all cells.
        for (int ix=0; ix<N; ++ix) {
            const int& current_cell=cptr[ix];
            if (current_cell >= ncells || current_cell < 0)  {
                throw std::runtime_error("chosen indices out of range");
            }

            searcher->find_neighbors(current_cell, threshold, false);
            if (collected.size()==0) {
                // Check here, otherwise median calculations fail.
                throw std::runtime_error("cell failed to count itself");
            }
            
            // Storing the identities of the cells (compressed).
            pack_index_vector(sorted_ids, collected.begin(), collected.end());
            SET_VECTOR_ELT(cellids, ix, allocVector(INTSXP, sorted_ids.size()));
            std::copy(sorted_ids.begin(), sorted_ids.end(), INTEGER(VECTOR_ELT(cellids, ix)));
            numptr[ix]=collected.size();
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
return mkString(e.what());
}

