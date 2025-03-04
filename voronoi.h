#ifndef H_VORONOI
#define H_VORONOI

#include "jc_voronoi.h"
#include <getopt.h>
#include <stdio.h>

typedef struct data_           data;
typedef struct parameter_      parameter;

struct parameter_
{
    jcv_real Ao;
    jcv_real Po;
    jcv_real Ka;
    jcv_real Kp;
};

struct data_
{
    jcv_diagram* diagram;
    jcv_point* positions;
    jcv_point* velocities;
    jcv_point* forces;
    const jcv_site* sites;

    jcv_real L, dL;
    jcv_real gamma_rate;
    jcv_real amount_of_def;
    jcv_real dt;
    parameter parameter;
    int N;
    int M;
    int N_pbc;
    int i;

    FILE* file;
};


void constantInit(int argc, char *argv[], data* sys);

#endif // H_VORONOI