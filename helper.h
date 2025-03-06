#ifndef HELPER_H
#define HELPER_H

#include "jc_voronoi.h"
#include <stdio.h>
#include "voronoi.h"

jcv_point jcv_add(const jcv_point a,const  jcv_point b);
jcv_point jcv_sub(const jcv_point a, const jcv_point b);
jcv_real jcv_cross(const jcv_point a, const jcv_point b);
jcv_real jcv_dot(const jcv_point a, const jcv_point b);
jcv_point jcv_mul(const jcv_real a, const jcv_point b);
jcv_point jcv_vec_mul(const jcv_point a, const jcv_point b);

jcv_real jcv_perimeter(const jcv_site* a);
jcv_real jcv_area(const jcv_site* a);

jcv_real jcv_lenght(const jcv_point* a);
jcv_real jcv_lenght_sq(const jcv_point* a);


void populate_points(jcv_point* points, int N, jcv_real L);
void pbc(jcv_point* p, jcv_real L, jcv_real gamma);
void addBoundary(data* sys, int i);

double drand(double min, double max);
void gaussian(jcv_real* n1, jcv_real* n2);

#endif // HELPER_H