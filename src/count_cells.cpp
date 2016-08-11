#include "objects.h"

extern "C" {

SEXP count_cells(SEXP exprs, SEXP distance, SEXP nsamp, SEXP sample_id, SEXP centers, SEXP cluster_info, SEXP curcells) try {
    finder fx(exprs, centers, cluster_info);

    // Checking samples.
    if (!isInteger(nsamp)|| LENGTH(nsamp)!=1) { throw std::runtime_error("number of samples must be an integer scalar"); }
    const int nsamples=asInteger(nsamp);
    if (nsamples <= 0) { throw std::runtime_error("number of samples must be positive"); }

    if (!isInteger(sample_id) || LENGTH(sample_id)!=int((fx.searcher->exprs).ncol)) {
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }
    const int* sample_ids=INTEGER(sample_id);
    for (int i=0; i<LENGTH(sample_id); ++i) {
        if (sample_ids[i] < 0 || sample_ids[i] >= nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
    }

    // Checking distances and chosen cells.
    if (!isReal(distance)|| LENGTH(distance)!=1) { throw std::runtime_error("distance must be a double-precision scalar"); }
    const double threshold=asReal(distance);
    if (!isInteger(curcells)) {
        throw std::runtime_error("chosen indices should be an integer vector");
    }
    const int N=LENGTH(curcells);
    const int ncells=fx.searcher->get_ncells();
    const size_t& nmarkers=fx.searcher->get_nmarkers();
    if (nmarkers==0) {
        throw std::runtime_error("number of markers should be positive");
    }
    const int* cptr=INTEGER(curcells);
    
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        std::deque<int*> count_ptrs;
        SET_VECTOR_ELT(output, 0, allocMatrix(INTSXP, N, nsamples));
        SEXP current=VECTOR_ELT(output, 0);
        count_ptrs.push_back(INTEGER(current));
        std::fill(count_ptrs.back(), count_ptrs.back()+N, 0);
        for (int inner_si=1; inner_si<nsamples; ++inner_si) { 
            count_ptrs.push_back(count_ptrs[inner_si-1] + N);
            std::fill(count_ptrs.back(), count_ptrs.back()+N, 0);
        }

        std::deque<double*> coord_ptrs;
        SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, N, nmarkers));
        current=VECTOR_ELT(output, 1);
        coord_ptrs.push_back(REAL(current));
        for (size_t mi=1; mi<nmarkers; ++mi) { 
            coord_ptrs.push_back(coord_ptrs[mi-1] + N);
        }

        size_t icx, mi;
        int si;   
        std::deque<double> intensities;
        const double* marker_exprs;
        size_t midpoint;
        bool iseven; 
        std::deque<size_t>& collected=(fx.searcher->neighbors);

        // Running through all cells.
        for (int ix=0; ix<N; ++ix) {
            const int& current_cell=cptr[ix];
            if (current_cell >= ncells || current_cell < 0)  {
                throw std::runtime_error("chosen indices out of range");
            }

            fx.searcher->find_neighbors(current_cell, threshold, false);
            if (collected.size()==0) {
                // Check here, otherwise median calculations fail.
                throw std::runtime_error("cell failed to count itself");
            }
            for (icx=0; icx<collected.size(); ++icx) {
                ++(count_ptrs[sample_ids[collected[icx]]][0]);
            }

            // Setting the medians.
            intensities.resize(collected.size());
            midpoint=intensities.size()/2;
            iseven=(intensities.size()%2==0);

            for (mi=0; mi<nmarkers; ++mi) {
                marker_exprs=(fx.searcher->exprs).dptr + mi;
                for (icx=0; icx<collected.size(); ++icx) {
                    intensities[icx]=marker_exprs[nmarkers*collected[icx]];
                }
                
                std::sort(intensities.begin(), intensities.end());
                coord_ptrs[mi][0]=(iseven ?
                        (intensities[midpoint-1]+intensities[midpoint])/2 :
                        (intensities[midpoint]));
                ++(coord_ptrs[mi]);
            }

            // Bumping forward the ptrs.
            for (si=0; si<nsamples; ++si) {
                ++(count_ptrs[si]);
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

}

