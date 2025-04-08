#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "voronoi.h"

void backwardEulerStep(data* sys);
void rk4Step(data* sys);
void rkf45Step(data* sys);
void eulerStep(data* sys);
void fireStep(data* sys);
void conjugateGradientStep(data* sys);
jcv_real line_search(data* sys, jcv_point* gradient, jcv_point *pk, jcv_real alpha, jcv_real rho, jcv_real c);
jcv_real line_search_local_minima(data* sys, jcv_point* gradient, jcv_point *pk, jcv_real alpha, jcv_real rho, jcv_real c);
void computeLandscape(data* sys);
#endif



        