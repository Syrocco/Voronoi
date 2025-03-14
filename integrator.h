#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "voronoi.h"

void eulerStep(data* sys);
void fireStep(data* sys);
void conjugateGradientStep(data* sys);
jcv_real line_search(data* sys, jcv_point* gradient, jcv_point *pk, double alpha, double rho, double c);
#endif



        