#include "objects.h"

SEXP get_knn_distance(SEXP coords, SEXP centers, SEXP clust_info, SEXP nn) try {
    if (!isInteger(nn) || LENGTH(nn)!=1) { 
        throw std::runtime_error("number of neighbours must be an integer scalar");
    }
    const size_t NN=asInteger(nn);
    finder fx(coords, centers, clust_info);
    const size_t nhyper=fx.searcher -> get_ncells();

    SEXP output=PROTECT(allocVector(REALSXP, nhyper));
    try {
        double* optr=REAL(output);
        std::deque<double>& distances=fx.searcher->distances;
        
        for (size_t h=0; h<nhyper; ++h) {
            fx.searcher->find_nearest_neighbors(h, NN, true);
            optr[h]=distances.back();
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

SEXP compute_density (SEXP coords, SEXP centers, SEXP clust_info, SEXP radius) try {
    if (!isReal(radius) || LENGTH(radius)!=1) { 
        throw std::runtime_error("radius must be a double-precision scalar"); 
    }
    const double rad=asReal(radius);
    finder fx(coords, centers, clust_info);
    const size_t nhyper=fx.searcher -> get_ncells();

    SEXP output=PROTECT(allocVector(REALSXP, nhyper));
    try {
        double* optr=REAL(output);
        std::deque<double>& distances=fx.searcher->distances;

        double diffdist;
        size_t xi;
        for (size_t h=0; h<nhyper; ++h) {
            fx.searcher->find_neighbors(h, rad, true);
            
            double& curdensity=(optr[h]=0);
            for (xi=0; xi<distances.size(); ++xi) {
                diffdist = 1 - std::pow(distances[xi]/rad, 3);
                curdensity += diffdist * diffdist * diffdist; // tricube weights.
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

