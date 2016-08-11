#ifndef CYDAR_H
#define CYDAR_H

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"


#include <stdexcept>
#include <algorithm>
#include <map>
#include <deque>
#include <bitset>
#include <cmath>

extern "C" {

SEXP get_knn_distance(SEXP, SEXP, SEXP, SEXP);

SEXP compute_density(SEXP, SEXP, SEXP, SEXP);

SEXP count_cells(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP find_counted(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_nndist(SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif 
