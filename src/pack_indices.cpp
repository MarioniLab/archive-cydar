#include "cydar.h"
#include "packer.hpp"

void unpack_indices (std::deque<int>& output, const int* start, const int* end) {
    output.clear();
    while (start!=end) {
        const int& curval=(*start);
        if (curval > 0) {
            if (!output.empty() && curval < output.back()) {
                throw std::runtime_error("absolute values of compressed indices must always increase");
            }
            output.push_back(curval);
        } else if (curval < 0) {
            if (output.empty() || *(start-1) < 0 || output.back() > -curval) {
                throw std::runtime_error("inappropriate negative values in compressed index vector");
            }
            while (output.back()!=-curval) {
                output.push_back(output.back()+1);
            }
        } else {
            throw std::runtime_error("zero values in compressed index vector");
        }
        ++start;
    }
    return;
}

SEXP pack_indices(SEXP assignments, SEXP compact) try {
    if (!isNewList(assignments)) { 
        throw std::runtime_error("assignments should be a list");
    }
    const int ngrps=LENGTH(assignments);

    if (!isLogical(compact) || LENGTH(compact)!=1) { 
        throw std::runtime_error("'compact' should be a logical scalar");
    }
    const bool compress=asLogical(compact);

    SEXP output=PROTECT(allocVector(VECSXP, ngrps)), current;
    try { 
        std::deque<int> sorted_ids, temp;
        const int* iptr=NULL;
        int ix, ndex;

        for (int g=0; g<ngrps; ++g) {
            current=VECTOR_ELT(assignments, g);
            if (!isInteger(current)) { 
                throw std::runtime_error("assignment vectors should be integer");
            }
            iptr=INTEGER(current);
            ndex=LENGTH(current);
            
            if (compress) {
                temp.assign(iptr, iptr+ndex);
                for (ix=0; ix<ndex; ++ix) { --(temp[ix]); } // getting to zero-indexing.
                pack_index_vector(sorted_ids, temp.begin(), temp.end());
                SET_VECTOR_ELT(output, g, allocVector(INTSXP, sorted_ids.size()));
                std::copy(sorted_ids.begin(), sorted_ids.end(), INTEGER(VECTOR_ELT(output, g)));
            } else {
                unpack_indices(temp, iptr, iptr+ndex);
                SET_VECTOR_ELT(output, g, allocVector(INTSXP, temp.size()));
                std::copy(temp.begin(), temp.end(), INTEGER(VECTOR_ELT(output, g)));
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

SEXP renew_indices (SEXP ref, SEXP compressed) try {
    if (!isNewList(compressed)) { 
        throw std::runtime_error("compressed should be a list");
    }
    const int ngrps=LENGTH(compressed);

    if (!isInteger(ref)) { 
        throw std::runtime_error("reference indices should be integer");
    }
    const int nref=LENGTH(ref);
    const int* rptr=INTEGER(ref);
    --rptr; // Get to 1-indexed.

    SEXP output=PROTECT(allocVector(VECSXP, ngrps)), current;
    try {
        std::deque<int> tmp, out;
        const int* iptr;
        int ndex;
        size_t i;

        for (int g=0; g<ngrps; ++g) {
            current=VECTOR_ELT(compressed, g);
            if (!isInteger(current)) { 
                throw std::runtime_error("assignment vectors should be integer");
            }
            iptr=INTEGER(current);
            ndex=LENGTH(current);

            // Unpacking to 1-indexed form, then converting to reference indices (zero-indexed).
            unpack_indices(tmp, iptr, iptr+ndex);
            for (i=0; i<tmp.size(); ++i) {
                int& curi=tmp[i];
                if (curi <= 0 || curi > nref) {
                    throw std::runtime_error("specified index outside the range of the reference vector");
                }
                curi=rptr[curi]-1;
            }

            // Repacking and saving.
            pack_index_vector(out, tmp.begin(), tmp.end());
            SET_VECTOR_ELT(output, g, allocVector(INTSXP, out.size()));
            std::copy(out.begin(), out.end(), INTEGER(VECTOR_ELT(output, g)));
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

