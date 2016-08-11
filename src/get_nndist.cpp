#include "objects.h"

extern "C" {

SEXP get_nndist(SEXP cells, SEXP centers, SEXP clust_info, SEXP nn, SEXP freq) try {
    if (!isInteger(nn) || LENGTH(nn)!=1) { 
        throw std::runtime_error("number of neighbours must be an integer scalar");
    }
    const size_t NN=asInteger(nn);
    if (!isInteger(freq) || LENGTH(freq)!=1) {
        throw std::runtime_error("downsampling frequency must be an integer scalar");
    }
    const size_t downsample=asInteger(freq);
    finder fx(cells, centers, clust_info);
    const size_t ncells=fx.searcher -> get_ncells();

    const size_t true_nrows=(ncells ? 1+int((ncells-1)/downsample) : 0);
    SEXP output=PROTECT(allocMatrix(REALSXP, NN, true_nrows));
    try {
        double* optr=REAL(output);
        std::deque<double>& distances=fx.searcher->distances;
        
        for (size_t h=0; h<ncells; h+=downsample) {
            fx.searcher->find_nearest_neighbors(h, NN+1, true); // Avoiding itself.
            std::copy(distances.begin()+1, distances.end(), optr);
            optr+=NN;
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
