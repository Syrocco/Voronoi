#ifndef H_INITIAL
#define H_INITIAL
#include "voronoi.h"

void randomInitial(data* sys);
void rsaInitial(data* sys, jcv_real min_distance);
void distribute_area(data* sys);

#endif