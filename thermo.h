#ifndef H_THERMO
#define H_THERMO

#include "voronoi.h"
#include "jc_voronoi.h"

void shear(data* sys);

jcv_real energy_total(const jcv_site* sites, const int N, const int N_pbc, const parameter* param);
jcv_real energy_unique(const jcv_site* site, const parameter* param);

void stress_total(data* sys, jcv_real stress[2][2]);
void stress_unique(const jcv_site* site, const parameter* param, jcv_real stress[2][2]);

void distance_moved(data* sys, jcv_point* old_positions, jcv_real* dist);

#endif