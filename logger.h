#ifndef H_LOGGER
#define H_LOGGER

#include <stdio.h>
#include "voronoi.h"

void loggers(data* sys);

void saveTXT(data* sys);
void computeThermo(data* sys);
void write(FILE* file, const char* filename, const jcv_point* points, const jcv_site* sites, int N, int N_pbc);


#endif