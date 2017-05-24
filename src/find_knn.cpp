#include "objects.h"

SEXP find_knn(SEXP cells, SEXP clust_centers, SEXP clust_info, SEXP nn, SEXP mode, SEXP query) try {
    if (!isInteger(nn) || LENGTH(nn)!=1) { 
        throw std::runtime_error("number of neighbours must be an integer scalar");
    }
    const size_t NN=asInteger(nn);
    if (NN<1) { 
        throw std::runtime_error("number of nearest neighbors must be positive");
    }
    auto searcher=generate_holder(cells, clust_centers, clust_info);
    const size_t nmarkers=searcher->get_nmarkers();

    // Just iterating across itself, if there are no query cells.
    const bool self_cycle=(query==R_NilValue);
    const double* qptr=NULL;
    size_t ncells;
    if (self_cycle) { 
        query=cells; 
        ncells=searcher->get_ncells();
    } else {
        matrix_info qcells=check_matrix(query);
        if (qcells.dptr==NULL) {
            throw std::runtime_error("host intensity matrix should be double-precision");
        }       
        qptr=qcells.dptr;
        if (qcells.nrow!=nmarkers) {
            throw std::runtime_error("host and target intensity matrix do not have same number of markers");
        }
        ncells=qcells.ncol;
    }

    // Getting the output mode.
    if (!isInteger(mode) || LENGTH(mode)!=1) { 
        throw std::runtime_error("mode should be an integer scalar");
    }
    int M=asInteger(mode);
    const bool keep_last=(M < 0);
    if (keep_last) { M*=-1; }
    const bool store_distances=(M%2==1), store_neighbors=(M>=2);

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        double* odptr=NULL;
        int *onptr=NULL;
        if (store_distances) { // Record all distances.
            if (keep_last) { 
                SET_VECTOR_ELT(output, 1, allocVector(REALSXP, ncells));
            } else {
                SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, NN, ncells));
            }
            odptr=REAL(VECTOR_ELT(output, 1));
        } else { 
            SET_VECTOR_ELT(output, 1, R_NilValue);
        }
        if (store_neighbors) { // Record all neighbours.
            if (keep_last) { 
                SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ncells));
            } else {
                SET_VECTOR_ELT(output, 0, allocMatrix(INTSXP, NN, ncells));
            }
            onptr=INTEGER(VECTOR_ELT(output, 0));
        } else { 
            SET_VECTOR_ELT(output, 0, R_NilValue);
        }
        
        // Iterating across cells, finding NNs and storing (last) distances or neighbors.
        std::deque<double>& distances=searcher->distances;
        std::deque<size_t>& neighbors=searcher->neighbors;
        std::deque<size_t>::const_iterator nIt;
        
        for (size_t h=0; h<ncells; ++h) {
            if (self_cycle) { 
                searcher->find_nearest_neighbors(h, NN, store_distances); 
            } else {
                searcher->find_nearest_neighbors(qptr, NN, store_distances); 
                qptr+=nmarkers;
            }

            if (store_distances) {
                if (keep_last) {
                    (*odptr)=distances.back();
                    ++odptr;
                } else {
                    std::copy(distances.begin(), distances.end(), odptr);
                    odptr+=NN;
                }
            }
            if (store_neighbors) {
                if (keep_last) { 
                    (*onptr)=neighbors.back();
                    ++onptr;
                } else {
                    for (nIt=neighbors.begin(); nIt!=neighbors.end(); ++nIt, ++onptr) {
                        (*onptr)=int(*nIt);
                    }
                }
            }
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);

    if (store_distances && store_neighbors) { 
        return output;
    } else if (store_distances) {
        return VECTOR_ELT(output, 1);
    } else if (store_neighbors) {
        return VECTOR_ELT(output, 0);
    } else {
        return R_NilValue;
    }
} catch (std::exception& e) {
    return mkString(e.what());
}


