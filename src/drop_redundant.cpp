#include "objects.h"

SEXP drop_redundant (SEXP actual_order, SEXP coords, SEXP centers, SEXP clust_info, SEXP threshold) try {
    auto searcher=generate_holder(coords, centers, clust_info);
    const size_t nhyper=searcher->get_ncells();
    const size_t nmarkers=searcher->get_nmarkers();

    if (!isReal(threshold) || LENGTH(threshold)!=1) { 
        throw std::runtime_error("threshold must be a double-precision scalar"); 
    }
    const double thresh=asReal(threshold);
    const double radius=thresh * std::sqrt(nmarkers);

    if (!isInteger(actual_order)) {
        throw std::runtime_error("actual_order order vector must be integer");
    } 
    if (LENGTH(actual_order)!=nhyper) {
        throw std::runtime_error("length of actual_order order vector must be equal to number of hyperspheres");
    }
    const int* ordering=INTEGER(actual_order);
    
    SEXP output=PROTECT(allocVector(LGLSXP, nhyper));
    try {
        int* optr=LOGICAL(output);
        std::fill(optr, optr+nhyper, 0);
        std::deque<size_t>& neighbors=searcher->neighbors;
        std::deque<bool> already_seen(nhyper, false);

        const double* dptr=(searcher->exprs).dptr;
        if (dptr==NULL) {
            throw std::runtime_error("intensities must be double-precision");
        }

        for (size_t h=0; h<nhyper; ++h) {
            const int& actual_index=ordering[h];
            if (already_seen[actual_index]) { continue; }
            optr[actual_index]=1;
            searcher->find_neighbors(actual_index, radius, false);
            const double* dptr_current=dptr + actual_index*nmarkers;

            for (size_t ni=0; ni<neighbors.size(); ++ni) {
                bool okay=false;
                const double* dptr_target=dptr + neighbors[ni] * nmarkers;
                for (size_t mi=0; mi<nmarkers; ++mi) {
                    if (std::abs(dptr_target[mi] - dptr_current[mi]) > thresh) {
                        okay=true;
                        break;
                    }
                }
                if (!okay) { 
                    already_seen[neighbors[ni]] = true;
                }
            }
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

