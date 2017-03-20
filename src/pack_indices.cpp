#include "cydar.h"
#include "packer.h"

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
                unpack_index_vector(temp, iptr, iptr+ndex);
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

