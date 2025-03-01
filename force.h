#ifndef FORCE_H
#define FORCE_H

#include "jc_voronoi.h"
jcv_point force_h(const jcv_real A, const jcv_real P, const jcv_point* h7, const jcv_point* h2, const jcv_point* h3);

jcv_point get_edge_force_ji(const jcv_site* si, const jcv_graphedge* edgei);
jcv_point get_edge_force_ii(const jcv_site* si);

jcv_real energy(const jcv_site* sites, const int N);
#endif
