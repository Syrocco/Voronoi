#ifndef FORCE_H
#define FORCE_H

#include "jc_voronoi.h"
#include "voronoi.h"

jcv_point force_h(const jcv_real A, const jcv_real P, const jcv_point* h7, const jcv_point* h2, const jcv_point* h3, const parameter* param);

jcv_point get_edge_force_ji(const jcv_site* si, const jcv_graphedge* edgei, const parameter* param);
jcv_point get_edge_force_ii(const jcv_site* si, const parameter* param);


void derivative(const jcv_point* ri, const jcv_point* rj, const jcv_point* rk, jcv_point jacobian[2]);

void compute_force(data* sys);

void shear(data* sys);

jcv_real energy_total(const jcv_site* sites, const int N, const int N_pbc, const parameter* param);
jcv_real energy_unique(const jcv_site* site, const parameter* param);

void stress_total(data* sys, jcv_real stress[2][2]);
void stress_unique(const jcv_site* site, const parameter* param, jcv_real stress[2][2]);

void distance_moved(data* sys, jcv_point* old_positions, jcv_real* dist);

#endif
