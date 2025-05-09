#ifndef H_VORONOI
#define H_VORONOI

#include "jc_voronoi.h"
#include <getopt.h>
#include <stdio.h>

typedef struct data_                 data;
typedef struct parameter_            parameter;
typedef struct logger_thermo_        logger_thermo;
typedef struct logger_snapshot_      logger_snapshot;
typedef struct memory_               memory;

#define NOISE 0

struct logger_thermo_{
    int n_log;
    int n_start;
    FILE* file;

    int compute_stress;
    int compute_dist_travelled;
};

struct logger_snapshot_{
    int n_log;
    int n_start;

    int include_boundary;
    FILE* file;

    int compute_stress;
    int compute_dist_travelled;
    int compute_area;
    int compute_perimeter;
};

struct memory_{
    jcv_point* old_positions_thermo;
    jcv_point* old_positions_snapshot;
    jcv_point* old_positions_strobo;
    int new;
};

struct parameter_
{
    jcv_real qo;
    jcv_real Ka;
    jcv_real Kp;
    jcv_real T;
    jcv_real gamma_rate;
    jcv_real rad;
    jcv_real k;
    
};

struct data_
{
    jcv_diagram* diagram;
    jcv_point* positions;
    jcv_point* velocities;
    jcv_point* forces;
    jcv_real* prefered_area;

    const jcv_site* sites;
    
    jcv_real L, dL; //dL = fraction of L to include at the boundaries for pbc
    
    jcv_real gamma;
    jcv_real gamma_max;
    int shear_start;
    int n_to_shear_max;
    int shear_cycle;
    
    jcv_real dt, dt_fire;
    parameter parameter;
    int N;
    int M;
    int N_pbc; //N + number of boundary points
    int i;
    jcv_real time;

    jcv_real n_frac_small;
    jcv_real size_large_over_small;
    
    logger_thermo info_thermo;
    logger_snapshot info_snapshot;
    logger_thermo info_strobo;

    memory old_info;
};


void constantInit(int argc, char *argv[], data* sys);

#endif // H_VORONOI