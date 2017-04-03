#ifndef OBJECTS_H
#define OBJECTS_H
#include "cydar.h"
#include "utils.h"
#include <queue>

struct naive_holder {
    naive_holder(SEXP, SEXP);
    virtual ~naive_holder();
   
    void find_neighbors(size_t, double, const bool);
    void find_neighbors(const double*, double, const bool);
    void find_nearest_neighbors(size_t, size_t, const bool);
    void find_nearest_neighbors(const double*, size_t, const bool);
    
    size_t get_ncells() const;
    size_t get_nmarkers() const;
    const std::deque<size_t>& get_used_markers() const;

    double compute_marker_distance(const double*, const double*) const;    
    virtual void search_all(const double*, double, const bool);
    virtual void search_nn (const double*, size_t, const bool);
  
    matrix_info exprs;
    std::deque<size_t> rows_to_use;
    std::deque<size_t> neighbors;
    std::deque<double> distances;
    typedef std::priority_queue<std::pair<double, int> > nearest;
    nearest current_nearest;
};

struct convex_holder : public naive_holder {
    convex_holder(SEXP, SEXP, SEXP, SEXP);
    ~convex_holder();
    void search_all(const double*, double, const bool);
    void search_nn (const double*, size_t, const bool);

    matrix_info centers;
    std::deque<int> clust_start;
    std::deque<int> clust_ncells;
    std::deque<const double*> clust_dist;
};

struct finder {
    finder (SEXP, SEXP, SEXP, SEXP);
    ~finder ();
    naive_holder* searcher;
};

#endif
