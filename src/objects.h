#ifndef OBJECTS_H
#define OBJECTS_H

#include "cydar.h"
#include "utils.h"
#include <queue>
#include <memory>

struct naive_holder {
    naive_holder(SEXP);
    virtual ~naive_holder();
   
    void find_neighbors(size_t, double, const bool);
    void find_neighbors(const double*, double, const bool);
    void find_nearest_neighbors(size_t, size_t, const bool);
    void find_nearest_neighbors(const double*, size_t, const bool);
    
    size_t get_ncells() const;
    size_t get_nmarkers() const;

    double compute_marker_sqdist(const double*, const double*) const;    
    virtual void search_all(const double*, double, const bool);
    virtual void search_nn (const double*, size_t, const bool);
  
    matrix_info exprs;
    std::deque<size_t> neighbors;
    std::deque<double> distances;
    typedef std::priority_queue<std::pair<double, int> > nearest;
    nearest current_nearest;
};

struct convex_holder : public naive_holder {
    convex_holder(SEXP, SEXP, SEXP);
    ~convex_holder();
    void search_all(const double*, double, const bool);
    void search_nn (const double*, size_t, const bool);

    matrix_info centers;
    std::deque<int> clust_start;
    std::deque<int> clust_ncells;
    std::deque<const double*> clust_dist;
};

std::unique_ptr<naive_holder> generate_holder(SEXP, SEXP, SEXP); 

#endif
