#ifndef H_THERMO
#define H_THERMO

#include "voronoi.h"
#include "jc_voronoi.h"

void shear(data* sys);

jcv_real energy_total(data* sys);
jcv_real energy_unique(const jcv_site* site, const jcv_real A, const parameter* param);
jcv_real harmonic_energy(const jcv_point r, const jcv_real size, const jcv_real k);

void stress_total(data* sys, jcv_real stress[2][2]);
void stress_unique_cell(data* sys, int num, jcv_real stress[2][2]);
void stress_unique_force(data* sys, int num, jcv_real stress[2][2]);

void distance_moved(data* sys, jcv_point* old_positions, jcv_real* dist);

#endif