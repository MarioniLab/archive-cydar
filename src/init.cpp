#include "cydar.h"
#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(compute_density, 4),
    REGISTER(compute_hyperstats, 4),
    REGISTER(count_cells, 5),
    REGISTER(get_nndist, 5),
    REGISTER(find_knn, 6),
    REGISTER(drop_redundant, 5),
    REGISTER(pack_indices, 2),
    REGISTER(recount_cells, 4),
    REGISTER(compute_median_int, 4),
    {NULL, NULL, 0}
};

void attribute_visible R_init_cydar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

