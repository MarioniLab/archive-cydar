#include "packer.hpp"
#include "cydar.h"
#include "objects.h"

/* For each original hypersphere, this function identifies the cells within the corresponding nested 
 * hypersphere; then identifies other nested hyperspheres which contain a subset of cells within the 
 * original hypersphere.
 */

SEXP recount_cells(SEXP exprs, SEXP use_markers, SEXP distance, SEXP centers, SEXP assignments, SEXP filter) try {
    finder fx(exprs, use_markers, R_NilValue, R_NilValue);
    const matrix_info& EXPRS=fx.searcher->exprs;
    const size_t& nmarkers=EXPRS.nrow;
    const size_t& ncells=EXPRS.ncol;
    const double* eptr=EXPRS.dptr;
    const std::deque<size_t>& used_markers=fx.searcher->get_used_markers();
    const size_t& nused=used_markers.size();

    // Get centers and assignments for each hypersphere.
    if (!isInteger(centers)) {
        throw std::runtime_error("'centers' must be an integer vector");
    }
    const int* cptr=INTEGER(centers);
    const int ncenters=LENGTH(centers);

    if (!isNewList(assignments)) { 
        throw std::runtime_error("'assignments' must be a list of packed assignments");
    } else if (LENGTH(assignments)!=ncenters) {
        throw std::runtime_error("length of 'assignments' must be equal to length of 'centers'");
    }
    
    if (!isReal(distance)|| LENGTH(distance)!=1) { throw std::runtime_error("distance must be a double-precision scalar"); }
    const double threshold=asReal(distance);
    if (!isInteger(filter)|| LENGTH(filter)!=1) { throw std::runtime_error("filter must be an integer scalar"); }
    const size_t filt=asInteger(filter);

    // Going through all centers and finding all cells that belong in the new space.
    std::deque<std::deque<int> > new_assignments(ncenters);
    std::deque<int> is_center(ncells, -1);
    std::deque<double> maxdist(ncenters);
    
    SEXP curass;
    const int* iptr=NULL;
    int ndex;
    size_t ix;
    std::deque<int> temp;

    size_t u;
    double dist, tmpdist;
    const double* centerpoint, *curpoint;

    for (int c=0; c<ncenters; ++c) { 
        curass=VECTOR_ELT(assignments, c);
        if (!isInteger(curass)) { 
            throw std::runtime_error("assignment vectors should be integer");
        }
        iptr=INTEGER(curass);
        ndex=LENGTH(curass);
        unpack_index_vector(temp, iptr, iptr+ndex);

        std::deque<int>& retained=new_assignments[c];
        centerpoint=eptr+nmarkers*(cptr[c]);
        for (ix=0; ix<temp.size(); ++ix) {
            int& curt=temp[ix];
            --curt; // become zero-indexed.

            curpoint=eptr+nmarkers*curt;
            dist=0;
            for (u=0; u<nused; ++u) {
                const size_t& m=used_markers[u];
                tmpdist=centerpoint[m] - curpoint[m];
                dist+=tmpdist*tmpdist;
            }
            dist=std::sqrt(dist);

            if (dist <= threshold) {
                retained.push_back(curt);
                if (dist > maxdist[c]) { maxdist[c]=dist; }
            }
        }

        if (retained.size() < filt) { retained.clear(); } // No point keeping it at all if it's too small.
        is_center[cptr[c]]=c; 
    }

    // Setting up output objects.
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, ncenters));
        SEXP cell_assignments=VECTOR_ELT(output, 0);
        SET_VECTOR_ELT(output, 1, allocVector(VECSXP, ncenters));
        SEXP cell_centers=VECTOR_ELT(output, 1);
        int* ccptr=NULL;

        int zerodex;
        size_t jx, kx;
        std::deque<bool> in_current_set(ncells, false);
        std::deque<std::deque<int> > all_subsets;
        std::deque<int> all_centers;
    
        SEXP outass;
        bool include_thyself;
        std::deque<int> packed;

        // Checking all cells within each hypersphere for whether they're centers of other hyperspheres.
        for (int c=0; c<ncenters; ++c) {
            curass=VECTOR_ELT(assignments, c);
            iptr=INTEGER(curass);
            ndex=LENGTH(curass);
            unpack_index_vector(temp, iptr, iptr+ndex);
            if (temp.size() < filt) {
                SET_VECTOR_ELT(cell_assignments, c, allocVector(VECSXP, 0));
                SET_VECTOR_ELT(cell_centers, c, allocVector(INTSXP, 0));
                continue;
            }

            for (kx=0; kx<temp.size(); ++kx) { 
                if (c==0) {Rprintf("%i\n", temp[kx]);  }
                in_current_set[--(temp[kx])]=true; // 1-indexing and listing those that are centers.
            }

            all_subsets.resize(1);
            for (ix=0; ix<temp.size(); ++ix) {
                const int& zerodex=temp[ix];
                const int& centredex=is_center[zerodex];
                if (centredex==-1 || centredex==c) { // Skipping if not center, or if current center. 
                    continue; 
                }  

                // Identifying the subset elements.
                const std::deque<int>& other_retained=new_assignments[centredex];
                if (other_retained.size() < filt) { continue; } // skipping if too small.
                std::deque<int>& host=all_subsets.back();
                kx=std::upper_bound(other_retained.begin(), other_retained.end(), temp.back()) - other_retained.begin();
                if (kx < other_retained.size()) { kx=other_retained.size(); }

                for (jx=0; jx<kx; ++jx) {
                    const int& curj=other_retained[jx];
                    if (in_current_set[curj]) { host.push_back(curj); }
                }
                
                // Discarding it, if it's not large enough.
                if (host.size() < filt) { 
                    host.clear();
                } else {
                    all_subsets.resize(all_subsets.size()+1);
                    all_centers.push_back(zerodex+1);
                }
            }

            all_subsets.pop_back();
            for (kx=0; kx<temp.size(); ++kx) { in_current_set[temp[kx]]=false; } // resetting for next round.
            std::deque<int>& retained=new_assignments[c];
            include_thyself=retained.size() >= filt;

            // Packing the output.
            SET_VECTOR_ELT(cell_assignments, c, allocVector(VECSXP, all_subsets.size()+include_thyself));
            outass=VECTOR_ELT(cell_assignments, c);
            SET_VECTOR_ELT(cell_centers, c, allocVector(INTSXP, all_subsets.size()+include_thyself));
            ccptr=INTEGER(VECTOR_ELT(cell_centers, c));
            std::copy(all_centers.begin(), all_centers.end(), ccptr);

            for (ix=0; ix<all_subsets.size(); ++ix) {
                std::deque<int>& curdices=all_subsets[ix];
                pack_index_vector(packed, curdices.begin(), curdices.end());
                SET_VECTOR_ELT(outass, ix, allocVector(INTSXP, packed.size()));
                std::copy(packed.begin(), packed.end(), INTEGER(VECTOR_ELT(outass, ix)));
            }

            // Including the nested hypersphere on the current center cell.
            if (include_thyself) {
                pack_index_vector(packed, retained.begin(), retained.end());
                SET_VECTOR_ELT(outass, all_subsets.size(), allocVector(INTSXP, packed.size()));
                std::copy(packed.begin(), packed.end(), INTEGER(VECTOR_ELT(outass, all_subsets.size())));
                ccptr[all_centers.size()]=cptr[c]+1;
            }

            all_subsets.clear();
            all_centers.clear();
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
