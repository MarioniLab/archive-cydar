#include "objects.h"

SEXP drop_redundant (SEXP coords, SEXP centers, SEXP clust_info, SEXP threshold) try {
    finder fx(coords, centers, clust_info);
    const size_t nhyper=fx.searcher -> get_ncells();
    const size_t nmarkers=fx.searcher -> get_nmarkers();

    if (!isReal(threshold) || LENGTH(threshold)!=1) { 
        throw std::runtime_error("threshold must be a double-precision scalar"); 
    }
    const double thresh=asReal(threshold);
    const double radius=thresh * std::sqrt(nmarkers);

    SEXP output=PROTECT(allocVector(LGLSXP, nhyper));
    try {
        int* optr=LOGICAL(output);
        std::fill(optr, optr+nhyper, 0);
        std::deque<size_t>& neighbors=fx.searcher->neighbors;
        std::deque<bool> already_seen(nhyper);

        size_t ni, mi;
        bool okay;
        const double* dptr=(fx.searcher -> exprs).dptr, *dptr_current=NULL, *dptr_target=NULL;
        if (dptr==NULL) {
            throw std::runtime_error("intensities must be double-precision");
        }

        for (size_t h=0; h<nhyper; ++h) {
            if (already_seen[h]) { continue; }
            optr[h]=1;
            fx.searcher->find_neighbors(h, radius, false);
            dptr_current=dptr + h*nmarkers;

            for (ni=0; ni<neighbors.size(); ++ni) {
                okay=false;
                dptr_target=dptr + neighbors[ni] * nmarkers;
                for (mi=0; mi<nmarkers; ++mi) {
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

