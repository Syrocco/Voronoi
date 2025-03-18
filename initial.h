#ifndef H_INITIAL
#define H_INITIAL
#include "voronoi.h"

void randomInitial(data* sys);
void rsaInitial(data* sys, jcv_real min_distance);
void distribute_area(data* sys);
void random_area(data* sys, jcv_real min_area, jcv_real max_area);
void read_from_dump_initial(data* sys, char* filename, int frame);
#endif