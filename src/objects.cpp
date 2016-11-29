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
    search(exprs.dptr + exprs.nrow*cell, false, 0, threshold, dist);
    return;
}

void naive_holder::find_neighbors (const double* current, double threshold, const bool dist) {
    search(current, false, 0, threshold, dist);
    return;
}

void naive_holder::find_nearest_neighbors (size_t cell, size_t nn, const bool dist) {
    if (cell >= exprs.ncol) { throw std::runtime_error("cell index out of range"); }
    search(exprs.dptr + exprs.nrow*cell, true, nn, 0, dist);
    return;
}

void naive_holder::find_nearest_neighbors (const double* current, size_t nn, const bool dist) {
    search(current, true, nn, 0, dist);
    return;
}

void naive_holder::search (const double* current, const bool nnsearch, size_t nn, double threshold, const bool dist) {
    neighbors.clear();
    distances.clear();

    const size_t& nmarkers=exprs.nrow;
    const size_t& ncells=exprs.ncol;
    const size_t& nuse=rows_to_use.size();
    const double* other=exprs.dptr;

    double curdist2=0, tmpdist;
    const double maxdist2=threshold*threshold;
    size_t u=0;
    for (size_t c=0; c<ncells; ++c, other+=nmarkers) {

        curdist2=0;
        for (u=0; u<nuse; ++u) {
            const size_t& m=rows_to_use[u];
            tmpdist = current[m] - other[m];
            curdist2 += tmpdist*tmpdist;
        }

        if (!nnsearch) { 
            if (curdist2 <= maxdist2) {
                neighbors.push_back(c);
                if (dist) {
                    distances.push_back(std::sqrt(curdist2));
                }
            }
        } else {
            if (current_nearest.size() < nn || curdist2 < current_nearest.top().first) {
                current_nearest.push(std::make_pair(curdist2, c));
                if (current_nearest.size() > nn) { 
                    current_nearest.pop();
                } 
            }
        }
    }

    if (nnsearch) { 
        // Converts information to neighbors/distances. Also clears 'nearest'.
        pqueue2deque(current_nearest, neighbors, dist, distances);
        // Converting back from squared distances.
        for (size_t i=0; i<distances.size(); ++i) { distances[i]=std::sqrt(distances[i]); }
    }
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

void convex_holder::search (const double* current, const bool nnsearch, size_t nn, double threshold, const bool dist) {
    neighbors.clear();
    distances.clear();
    const size_t& nmarkers=exprs.nrow;
    const size_t& ncenters=centers.ncol;
    const double* centerx=centers.dptr;
    if (nnsearch) { threshold = R_PosInf; }
    const size_t& nuse=rows_to_use.size();

    // More temporaries.
    double dist2center, dist2cell, lower_bd, tmp;
    size_t u=0;
    int index;
    const double* other;

    // Computing the distance to each center, and deciding whether to proceed for each cluster.
    for (size_t center=0; center<ncenters; ++center, centerx+=nmarkers) {
        dist2center=0;
        for (u=0; u<nuse; ++u) {
            const size_t& m=rows_to_use[u];
            tmp=centerx[m] - current[m];
            dist2center+=tmp*tmp;
        }
        dist2center=std::sqrt(dist2center);

        const int& cur_ncells=clust_ncells[center];
        if (!cur_ncells) { continue; }
        const double* cur_dist=clust_dist[center];
        const double& maxdist=cur_dist[cur_ncells-1];

        if (threshold + maxdist > dist2center) {
            // Cells within this cluster are potentially countable; proceeding to count them
            lower_bd=dist2center - threshold;
//            upper_bd=dist2center + threshold; // Doesn't help much.
            index=std::lower_bound(cur_dist, cur_dist + cur_ncells, lower_bd)-cur_dist;
            const int& cur_start=clust_start[center];
            other=exprs.dptr + nmarkers * (cur_start + index);

            for (; index<cur_ncells; ++index, other+=nmarkers) {
                dist2cell=0;
                for (u=0; u<nuse; ++u) {
                    const size_t& m=rows_to_use[u];
                    tmp=current[m] - other[m];
                    dist2cell+=tmp*tmp;
                }
                dist2cell=std::sqrt(dist2cell);
                    
                if (!nnsearch) {
                    if (dist2cell <= threshold) {
                        neighbors.push_back(cur_start + index);
                        if (dist) {
                            distances.push_back(dist2cell);
                        }
                    } 
                } else {
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
        }
    }

    if (nnsearch) { 
        // Converts information to neighbors/distances. Also clears 'nearest'.
        pqueue2deque(current_nearest, neighbors, dist, distances);
    }
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

