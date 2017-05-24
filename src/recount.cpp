#include "packer.h"
#include "cydar.h"
#include "objects.h"

/* For each original hypersphere, this function identifies the cells within the corresponding nested 
 * hypersphere; then identifies other nested hyperspheres which contain a subset of cells within the 
 * original hypersphere.
 */

SEXP recount_cells(SEXP exprs, SEXP distance, SEXP centers, SEXP assignments) try {
    const matrix_info& EXPRS=check_matrix(exprs);
    const size_t& nmarkers=EXPRS.nrow;
    const double* eptr=EXPRS.dptr;

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
    const double threshold2=asReal(distance) * asReal(distance);

    // Setting up output objects.
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, ncenters));
        SEXP cell_assignments=VECTOR_ELT(output, 0);
        SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ncenters));
        int* countptr=INTEGER(VECTOR_ELT(output, 1));

        // Going through all centers and finding all cells that belong in the new space.
        std::deque<std::deque<int> > new_assignments(ncenters);
        std::deque<int> temp, packed;
        
        for (int c=0; c<ncenters; ++c) { 
            SEXP curass=VECTOR_ELT(assignments, c);
            if (!isInteger(curass)) { 
                throw std::runtime_error("assignment vectors should be integer");
            }
            const int* iptr=INTEGER(curass);
            const int ndex=LENGTH(curass);
            unpack_index_vector(temp, iptr, iptr+ndex);
            
            std::deque<int>& retained=new_assignments[c];
            const double* centerpoint=eptr+nmarkers*(cptr[c]);
            for (size_t ix=0; ix<temp.size(); ++ix) {
                int& curt=temp[ix];
                --curt; // become zero-indexed.
                
                const double* curpoint=eptr+nmarkers*curt;
                double dist2=0, tmpdist;
                for (size_t m=0; m<nmarkers; ++m) {
                    tmpdist=centerpoint[m] - curpoint[m];
                    dist2+=tmpdist*tmpdist;
                }
                if (dist2 <= threshold2) { retained.push_back(curt); }
            }

            // Storing the output.
            countptr[c]=retained.size();
            pack_index_vector(packed, retained.begin(), retained.end());
            SET_VECTOR_ELT(cell_assignments, c, allocVector(INTSXP, packed.size()));
            std::copy(packed.begin(), packed.end(), INTEGER(VECTOR_ELT(cell_assignments, c)));
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
