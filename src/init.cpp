#include "cydar.h"
#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(get_knn_distance, 4),
    REGISTER(compute_density, 4),
    REGISTER(count_cells, 7),
    REGISTER(find_counted, 5),
    REGISTER(get_nndist, 5),
    {NULL, NULL, 0}
};

void attribute_visible R_init_cydar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

