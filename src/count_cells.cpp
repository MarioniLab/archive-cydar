#include "objects.h"
#include "packer.hpp"

SEXP count_cells(SEXP exprs, SEXP distance, SEXP nsamp, SEXP sample_id, SEXP centers, SEXP cluster_info, SEXP curcells, SEXP curmarks) try {
    finder fx(exprs, curmarks, centers, cluster_info);

    // Checking samples and computing sample weights.
    if (!isInteger(nsamp)|| LENGTH(nsamp)!=1) { throw std::runtime_error("number of samples must be an integer scalar"); }
    const int nsamples=asInteger(nsamp);
    if (nsamples <= 0) { throw std::runtime_error("number of samples must be positive"); }

    if (!isInteger(sample_id) || LENGTH(sample_id)!=int((fx.searcher->exprs).ncol)) {
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }
    const int* sample_ids=INTEGER(sample_id);
    double* sample_weights=(double*)R_alloc(nsamples, sizeof(double));
    std::fill(sample_weights, sample_weights+nsamples, 0);
    for (int i=0; i<LENGTH(sample_id); ++i) {
        const int& cursample=sample_ids[i];
        if (cursample < 0 || cursample >= nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
        ++(sample_weights[cursample]);
    }
    for (int i=0; i<nsamples; ++i) { // Reciprocal of the total number of cells.
        sample_weights[i]=1/sample_weights[i];
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
    const std::deque<size_t>& used_markers=fx.searcher->get_used_markers();
    const size_t& nused=used_markers.size();
   
    // Setting up output vectors. 
    SEXP output=PROTECT(allocVector(VECSXP, 3));
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
        std::fill(coord_ptrs.back(), coord_ptrs.back()+N, R_NaReal);
        for (size_t mi=1; mi<nmarkers; ++mi) {
            coord_ptrs.push_back(coord_ptrs[mi-1]+N);
            std::fill(coord_ptrs.back(), coord_ptrs.back()+N, R_NaReal);
        }

        SET_VECTOR_ELT(output, 2, allocVector(VECSXP, N));
        SEXP cellids=VECTOR_ELT(output, 2);
        std::deque<int> sorted_ids;

        size_t icx, u;
        int si;   
        std::deque<std::pair<double, int> > intensities;
        const double* marker_exprs;
        size_t midpoint;
        double total_weight, midweight, cumweight;
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
            
            // Storing the identities of the cells (compressed).
            pack_index_vector(sorted_ids, collected.begin(), collected.end());
            SET_VECTOR_ELT(cellids, ix, allocVector(INTSXP, sorted_ids.size()));
            std::copy(sorted_ids.begin(), sorted_ids.end(), INTEGER(VECTOR_ELT(cellids, ix)));
            
            // Computing counts and total weights.
            total_weight=0;
            for (icx=0; icx<collected.size(); ++icx) {
                const int& cursample=sample_ids[collected[icx]];
                ++(count_ptrs[cursample][0]);
                total_weight+=sample_weights[cursample];
            }

            // Setting the weighted medians (to avoid large samples from dominating the location).
            intensities.resize(collected.size());
            midweight=total_weight/2;
            for (u=0; u<nused; ++u) {
                const size_t& mi=used_markers[u];

                marker_exprs=(fx.searcher->exprs).dptr + mi;
                for (icx=0; icx<collected.size(); ++icx) {
                    const int& curneighbor=collected[icx];
                    intensities[icx].first=marker_exprs[nmarkers*curneighbor];
                    intensities[icx].second=sample_ids[curneighbor];
                }
                
                std::sort(intensities.begin(), intensities.end());
                cumweight=0;
                for (midpoint=0; midpoint<intensities.size(); ++midpoint) {
                    cumweight += sample_weights[intensities[midpoint].second];
                    if (cumweight/total_weight >= 0.5) { break; }
                }
               
                if (midpoint==intensities.size()) {
                    // Only possible if total_weights is zero. 
                    coord_ptrs[mi][0]=R_NaReal;
                } else {
                    if (cumweight/total_weight==0.5) {
                        coord_ptrs[mi][0]=(intensities[midpoint].first + intensities[midpoint+1].first)/2;
                    } else {
                        coord_ptrs[mi][0]=intensities[midpoint].first;                
                    }
                }
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

