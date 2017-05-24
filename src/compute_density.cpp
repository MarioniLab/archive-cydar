#include "objects.h"

SEXP compute_density (SEXP coords, SEXP centers, SEXP clust_info, SEXP radius) try {
    if (!isReal(radius) || LENGTH(radius)!=1) { 
        throw std::runtime_error("radius must be a double-precision scalar"); 
    }
    const double rad=asReal(radius);
    auto searcher=generate_holder(coords, centers, clust_info);
    const size_t nhyper=searcher->get_ncells();

    SEXP output=PROTECT(allocVector(REALSXP, nhyper));
    try {
        double* optr=REAL(output);
        std::deque<double>& distances=searcher->distances;

        for (size_t h=0; h<nhyper; ++h) {
            searcher->find_neighbors(h, rad, true);
            
            double& curdensity=(optr[h]=0);
            for (size_t xi=0; xi<distances.size(); ++xi) {
                double diffdist = 1 - std::pow(distances[xi]/rad, 3);
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

