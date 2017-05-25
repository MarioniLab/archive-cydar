#include "objects.h"
    
typedef std::priority_queue<std::pair<double, int> > nearest;

void pqueue2deque(nearest& collected_cells, std::deque<size_t>& final_cells, const bool dist, std::deque<double>& distances) {
    while (!collected_cells.empty()) {
        final_cells.push_front(collected_cells.top().second);
        if (dist) {
            distances.push_front(std::sqrt(collected_cells.top().first)); // rooting, as distances stored as squares.
        }
        collected_cells.pop();
    }
    return;
}

/****************** Naive search object *********************/

naive_holder::naive_holder (SEXP ex) : exprs(check_matrix(ex)) {}

naive_holder::~naive_holder() { }

size_t naive_holder::get_ncells() const { return exprs.ncol; }

size_t naive_holder::get_nmarkers() const { return exprs.nrow; }

void naive_holder::find_neighbors (size_t cell, double threshold, const bool dist) {
    if (cell >= exprs.ncol) { throw std::runtime_error("cell index out of range"); }
    search_all(exprs.dptr + exprs.nrow*cell, threshold, dist);
    return;
}

void naive_holder::find_neighbors (const double* current, double threshold, const bool dist) {
    search_all(current, threshold, dist);
    return;
}

void naive_holder::find_nearest_neighbors (size_t cell, size_t nn, const bool dist) {
    if (cell >= exprs.ncol) { throw std::runtime_error("cell index out of range"); }
    search_nn(exprs.dptr + exprs.nrow*cell, nn+1, dist);

    // Removing the cell itself, if it's in the NN range. Otherwise removing the last NN.
    size_t i=0;
    for (; i<neighbors.size()-1; ++i) {
        if (neighbors[i]==cell) { break; }
    }
    neighbors.erase(neighbors.begin()+i);
    if (dist) {
        distances.erase(distances.begin()+i);
    }
    return;
}

void naive_holder::find_nearest_neighbors (const double* current, size_t nn, const bool dist) {
    search_nn(current, nn, dist);
    return;
}

double naive_holder::compute_marker_sqdist(const double* x, const double* y) const {
    double out=0;
    for (size_t m=0; m<exprs.nrow; ++m) {
        double tmp=x[m]-y[m];
        out+=tmp*tmp;
    }
    return out;
}

void naive_holder::search_all(const double* current, double threshold, const bool dist) {
    neighbors.clear();
    distances.clear();

    const size_t& nmarkers=exprs.nrow;
    const size_t& ncells=exprs.ncol;
    const double* other=exprs.dptr;
    const double threshold2=threshold*threshold; // squaring.

    double curdist2=0;
    for (size_t c=0; c<ncells; ++c, other+=nmarkers) {
        curdist2=compute_marker_sqdist(current, other);
        if (curdist2 <= threshold2) {
            neighbors.push_back(c);
            if (dist) {
                distances.push_back(std::sqrt(curdist2));
            }
        }
    }
    return;
} 

void naive_holder::search_nn (const double* current, size_t nn, const bool dist) {
    neighbors.clear();
    distances.clear();

    const size_t& nmarkers=exprs.nrow;
    const size_t& ncells=exprs.ncol;
    const double* other=exprs.dptr;

    double curdist2=0;
    for (size_t c=0; c<ncells; ++c, other+=nmarkers) {
        curdist2=compute_marker_sqdist(current, other);
        if (current_nearest.size() < nn || curdist2 < current_nearest.top().first) {
            current_nearest.push(std::make_pair(curdist2, c));
            if (current_nearest.size() > nn) { 
                current_nearest.pop();
            } 
        }
    }

    // Converts information to neighbors/distances. Also clears 'nearest'.
    pqueue2deque(current_nearest, neighbors, dist, distances);
    return;
} 

/****************** Convex search object *********************/

convex_holder::convex_holder(SEXP ex, SEXP cen, SEXP info) : naive_holder(ex), centers(check_matrix(cen)) {
    const size_t& ncenters=centers.ncol;
    for (size_t i=0; i<ncenters; ++i) {
        SEXP current=VECTOR_ELT(info, i);
        if (!isNewList(current) || LENGTH(current)!=2) {
            throw std::runtime_error("list elements must be of length 2");
        }

        SEXP start=VECTOR_ELT(current, 0);
        if (!isInteger(start) || LENGTH(start)!=1) { 
            throw std::runtime_error("starting ID must be an integer scalar");
        }
        clust_start.push_back(asInteger(start));

        SEXP distances=VECTOR_ELT(current, 1);
        if (!isReal(distances)) { 
            throw std::runtime_error("distances must be a double-precision vector");
        } 
        clust_dist.push_back(REAL(distances));
        clust_ncells.push_back(LENGTH(distances));
    }
    return;
}

convex_holder::~convex_holder() { }

void convex_holder::search_all (const double* current, double threshold, const bool dist) {
    neighbors.clear();
    distances.clear();
    const size_t& nmarkers=exprs.nrow;
    const size_t& ncenters=centers.ncol;
    const double* centerx=centers.dptr;
    const double threshold2=threshold*threshold; // squaring.

    // Computing the distance to each center, and deciding whether to proceed for each cluster.
    for (size_t center=0; center<ncenters; ++center, centerx+=nmarkers) {
        const int& cur_ncells=clust_ncells[center];
        if (!cur_ncells) { continue; }

        const double dist2center=std::sqrt(compute_marker_sqdist(current, centerx));
        const double* cur_dist=clust_dist[center];
        const double& maxdist=cur_dist[cur_ncells-1];
        if (threshold + maxdist < dist2center) { continue; }

        /* Cells within this cluster are potentially countable; jumping to the first countable cell,
         * according to the triangle inequality. We could also define the last countable cell by the 
         * reverse triangle inequality, but the clusters are too compact for that to come into play.
         */
        const double lower_bd=dist2center-threshold;
        const int firstcell=std::lower_bound(cur_dist, cur_dist + cur_ncells, lower_bd) - cur_dist;
//        const double upper_bd=dist2center + threshold;
//        const int lastcell=std::upper_bound(cur_dist + firstcell, cur_dist + cur_ncells, upper_bd) - cur_dist;
        
        const int& cur_start=clust_start[center];
        const double* other=exprs.dptr + nmarkers * (cur_start + firstcell);
        for (int index=firstcell; index<cur_ncells; ++index, other+=nmarkers) {

            const double dist2cell2=compute_marker_sqdist(current, other);
            if (dist2cell2 <= threshold2) {
                neighbors.push_back(cur_start + index);
                if (dist) {
                    distances.push_back(std::sqrt(dist2cell2));
                }
            } 
        }
    }
    return;
}  

void convex_holder::search_nn(const double* current, size_t nn, const bool dist) {
    neighbors.clear();
    distances.clear();
    const size_t& nmarkers=exprs.nrow;
    const size_t& ncenters=centers.ncol;
    const double* centerx=centers.dptr;
    double threshold2 = R_PosInf;

    /* Computing distances to all centers and sorting them.
     * The aim is to go through the nearest centers first, to 
     * get the shortest 'threshold' possible.
     */
    std::deque<std::pair<double, size_t> > center_order(ncenters); 
    for (size_t center=0; center<ncenters; ++center, centerx+=nmarkers) {
        center_order[center].first=std::sqrt(compute_marker_sqdist(current, centerx));
        center_order[center].second=center;
    }
    std::sort(center_order.begin(), center_order.end());

    // Computing the distance to each center, and deciding whether to proceed for each cluster.
    for (size_t cx=0; cx<ncenters; ++cx) {
        const size_t& center=center_order[cx].second;
        const double& dist2center=center_order[cx].first;
        const double* centerx=centers.dptr+nmarkers*center;

        const int& cur_ncells=clust_ncells[center];
        if (!cur_ncells) { continue; }
        const double* cur_dist=clust_dist[center];
        const double& maxdist=cur_dist[cur_ncells-1];

        int firstcell=0;
//        double upper_bd=R_PosInf;
        if (R_FINITE(threshold2)) {
            const double threshold = std::sqrt(threshold2);
            if (threshold + maxdist < dist2center) { continue; }

            // Cells within this cluster are potentially countable; proceeding to count them
            const double lower_bd=dist2center - threshold;
            firstcell=std::lower_bound(cur_dist, cur_dist + cur_ncells, lower_bd)-cur_dist;
//            upper_bd = threshold + dist2center;
        }

        const int& cur_start=clust_start[center];
        const double* other=exprs.dptr + nmarkers * (cur_start + firstcell);
        for (int index=firstcell; index<cur_ncells; ++index, other+=nmarkers) {
//            if (cur_dist[index] > upper_bd) { 
//                break; 
//            }

            const double dist2cell2=compute_marker_sqdist(current, other);                   
            if (current_nearest.size() < nn || dist2cell2 <= threshold2) {
                current_nearest.push(std::make_pair(dist2cell2, cur_start + index));
                if (current_nearest.size() > nn) { 
                    current_nearest.pop();
                } 
                if (current_nearest.size()==nn) {
                    threshold2=current_nearest.top().first; // Shrinking the threshold, if an earlier NN has been found.
//                    upper_bd=std::sqrt(threshold2) + dist2center;
                } 
            }
        }
    }

    // Converts information to neighbors/distances. Also clears 'nearest'.
    pqueue2deque(current_nearest, neighbors, dist, distances);
    return;
}  

/****************** Finder *********************/

std::unique_ptr<naive_holder> generate_holder(SEXP coords, SEXP centers, SEXP clust_info) {
    if (centers==R_NilValue || clust_info==R_NilValue) { 
        return std::unique_ptr<naive_holder>(new naive_holder(coords));
    } else {
        return std::unique_ptr<naive_holder>(new convex_holder(coords, centers, clust_info));
    }
}

