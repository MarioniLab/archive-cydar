#include "objects.h"
    
typedef std::priority_queue<std::pair<double, int> > nearest;

void pqueue2deque(nearest& collected_cells, std::deque<size_t>& final_cells,
        const bool dist, std::deque<double>& distances) {
    while (!collected_cells.empty()) {
        final_cells.push_front(collected_cells.top().second);
        if (dist) {
            distances.push_front(collected_cells.top().first);
        }
        collected_cells.pop();
    }
    return;
}

/****************** Naive search object *********************/

naive_holder::naive_holder (SEXP ex, SEXP use) : exprs(check_matrix(ex)) {

    if (use!=R_NilValue) {
        if (!isLogical(use)) {
            throw std::runtime_error("'use' vector should be logical or NULL");
        }
        const int* uptr=LOGICAL(use);
        for (int u=0; u<LENGTH(use); ++u) {
            if (uptr[u]) { rows_to_use.push_back(u); }
        }
    } else {
        for (int u=0; u<exprs.nrow; ++u) {
            rows_to_use.push_back(u);
        }
    }
    return;
}

naive_holder::~naive_holder() { }

size_t naive_holder::get_ncells() const { return exprs.ncol; }

size_t naive_holder::get_nmarkers() const { return exprs.nrow; }

const std::deque<size_t>& naive_holder::get_used_markers() const { return rows_to_use; }

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

double naive_holder::compute_marker_distance(const double* x, const double* y) const {
    double tmp, out=0;
    const size_t& nuse=rows_to_use.size();
    for (size_t u=0; u<nuse; ++u) {
        const size_t& m=rows_to_use[u];
        tmp=x[m]-y[m];
        out+=tmp*tmp;
    }
    return std::sqrt(out);
}

void naive_holder::search_all(const double* current, double threshold, const bool dist) {
    neighbors.clear();
    distances.clear();

    const size_t& nmarkers=exprs.nrow;
    const size_t& ncells=exprs.ncol;
    const double* other=exprs.dptr;

    double curdist=0;
    for (size_t c=0; c<ncells; ++c, other+=nmarkers) {
        curdist=compute_marker_distance(current, other);
        if (curdist <= threshold) {
            neighbors.push_back(c);
            if (dist) {
                distances.push_back(curdist);
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

    double curdist=0;
    for (size_t c=0; c<ncells; ++c, other+=nmarkers) {
        curdist=compute_marker_distance(current, other);
        if (current_nearest.size() < nn || curdist < current_nearest.top().first) {
            current_nearest.push(std::make_pair(curdist, c));
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

convex_holder::convex_holder(SEXP ex, SEXP use, SEXP cen, SEXP info) : naive_holder(ex, use), centers(check_matrix(cen)) {
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

    // More temporaries.
    double dist2center, dist2cell, lower_bd;
    int index;
    const double* other;

    // Computing the distance to each center, and deciding whether to proceed for each cluster.
    for (size_t center=0; center<ncenters; ++center, centerx+=nmarkers) {
        const int& cur_ncells=clust_ncells[center];
        if (!cur_ncells) { continue; }

        dist2center=compute_marker_distance(current, centerx);
        const double* cur_dist=clust_dist[center];
        const double& maxdist=cur_dist[cur_ncells-1];
        if (threshold + maxdist < dist2center) { continue; }

        // Cells within this cluster are potentially countable; proceeding to count them
        lower_bd=dist2center - threshold;
//        upper_bd=dist2center + threshold; // Doesn't help much.
        index=std::lower_bound(cur_dist, cur_dist + cur_ncells, lower_bd)-cur_dist;
        const int& cur_start=clust_start[center];
        other=exprs.dptr + nmarkers * (cur_start + index);

        for (; index<cur_ncells; ++index, other+=nmarkers) {
            dist2cell=compute_marker_distance(current, other);
            if (dist2cell <= threshold) {
                neighbors.push_back(cur_start + index);
                if (dist) {
                    distances.push_back(dist2cell);
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
    double threshold = R_PosInf;

    /* Computing distances to all centers and sorting them.
     * The aim is to go through the nearest centers first, to 
     * get the shortest 'threshold' possible.
     */
    std::deque<std::pair<double, size_t> > center_order(ncenters); 
    size_t center;
    for (center=0; center<ncenters; ++center, centerx+=nmarkers) {
        center_order[center].first=compute_marker_distance(current, centerx);
        center_order[center].second=center;
    }
    std::sort(center_order.begin(), center_order.end());

    // More temporaries.
    double dist2center, dist2cell, lower_bd;
    int index;
    const double* other;

    // Computing the distance to each center, and deciding whether to proceed for each cluster.
    for (size_t cx=0; cx<ncenters; ++cx) {
        center=center_order[cx].second;
        dist2center=center_order[cx].first;
        centerx=centers.dptr+nmarkers*center;

        const int& cur_ncells=clust_ncells[center];
        if (!cur_ncells) { continue; }
        const double* cur_dist=clust_dist[center];
        const double& maxdist=cur_dist[cur_ncells-1];

        if (R_FINITE(threshold)) {
            if (threshold + maxdist < dist2center) { continue; }
            // Cells within this cluster are potentially countable; proceeding to count them
            lower_bd=dist2center - threshold;
//          upper_bd=dist2center + threshold; // Doesn't help much.
            index=std::lower_bound(cur_dist, cur_dist + cur_ncells, lower_bd)-cur_dist;
        } else {
            index=0;
        }
        const int& cur_start=clust_start[center];
        other=exprs.dptr + nmarkers * (cur_start + index);

        for (; index<cur_ncells; ++index, other+=nmarkers) {
            dist2cell=compute_marker_distance(current, other);
                    
            if (current_nearest.size() < nn || dist2cell <= threshold) {
                current_nearest.push(std::make_pair(dist2cell, cur_start + index));
                if (current_nearest.size() > nn) { 
                    current_nearest.pop();
                } 
                if (current_nearest.size()==nn) {
                    threshold=current_nearest.top().first; // Shrinking the threshold, if an earlier NN has been found.
                } 
            }
        }
    }

    // Converts information to neighbors/distances. Also clears 'nearest'.
    pqueue2deque(current_nearest, neighbors, dist, distances);
    return;
}  

/****************** Finder *********************/

finder::finder (SEXP coords, SEXP markers, SEXP centers, SEXP clust_info) {
    if (centers==R_NilValue || clust_info==R_NilValue) { 
        searcher=new naive_holder(coords, markers);
    } else {
        searcher=new convex_holder(coords, markers, centers, clust_info);
    }
    return;
}

finder::~finder() {
    delete searcher;
    return;
}

