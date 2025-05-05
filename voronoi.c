#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"
#include "mersenne.c"
#include "helper.h"
#include "force.h"
#include "initial.h"
#include "voronoi.h"
#include "integrator.h"
#include <getopt.h>
#include <errno.h>
#include "thermo.h"
#include "logger.h"

int main(int argc, char *argv[]){
    init_genrand(0);

    data sys;
    sys.parameter.qo = 3.5;
    sys.parameter.Ka = 0;
    sys.parameter.Kp = 1;
    sys.parameter.rad = 0.25;
    sys.parameter.k = 2;

    sys.size_large_over_small = 3.0/3.0;
    sys.n_frac_small = 0.5;

    sys.N = 300;
    sys.M = 200000000; 
    sys.dt = 0.001;
    sys.parameter.T = 0.0;
    sys.parameter.gamma_rate = 0.0001;
    sys.gamma_max = 3;
    sys.shear_start = 1;
    sys.dt_fire = 0.1;
    sys.shear_cycle = 0;

    sys.info_snapshot.n_log = 100;
    sys.info_snapshot.n_start = 0;
    sys.info_snapshot.include_boundary = 0;
    sys.info_snapshot.compute_stress = 1;
    sys.info_snapshot.compute_dist_travelled = 0;
    sys.info_snapshot.compute_area = 1;
    sys.info_snapshot.compute_perimeter = 1;

    sys.info_thermo.n_log = 1;
    sys.info_thermo.n_start = 0;
    sys.info_thermo.compute_stress = 1;
    sys.info_thermo.compute_dist_travelled = 0;

    sys.info_strobo.compute_stress = 1;
    sys.info_strobo.compute_dist_travelled = 1;

    
    constantInit(argc, argv, &sys);
    jcv_diagram diagram;
    memset(&(diagram), 0, sizeof(jcv_diagram));
    jcv_point positions[9*sys.N]; //9*sys.N == max N_pbc
    jcv_point velocities[sys.N];
    jcv_point forces[sys.N];
    jcv_real prefered_area[9*sys.N];
    sys.diagram = &diagram;
    sys.positions = positions;
    sys.velocities = velocities;
    sys.forces = forces;
    sys.prefered_area = prefered_area;
    

    jcv_point old_positions1[sys.N];
    sys.old_info.old_positions_thermo = old_positions1;

    jcv_point old_positions2[sys.N];
    sys.old_info.old_positions_snapshot = old_positions2;

    jcv_point old_positions3[sys.N];
    sys.old_info.old_positions_strobo = old_positions3;

    //randomInitial(&sys);
    distribute_area(&sys);
    //random_area(&sys, 0.5, 1.5);
    rsaInitial(&sys, 0.4);
    
    //conjugateGradientStep(&sys);
    fireStep(&sys);
    for (sys.i = 0; sys.i < sys.M; sys.i++){
        
        //backwardEulerStep(&sys);
        rkf45Step(&sys);
        //rk4Step(&sys);
        //eulerStep(&sys);
        //fireStep(&sys);
        //conjugateGradientStep(&sys);
        
    }
    jcv_diagram_generate(sys.N_pbc, sys.positions, NULL, NULL, sys.diagram);
    sys.sites = jcv_diagram_get_sites(sys.diagram);
    compute_force(&sys);
    computeThermo(&sys);
    saveTXT(&sys);
    //check_force(&sys);
    /*  jcv_diagram_generate(sys.N_pbc, sys.positions, NULL, NULL, sys.diagram);
    sys.sites = jcv_diagram_get_sites(sys.diagram);
    compute_force(&sys);
    computeThermo(&sys);
    computeLandscape(&sys);
    saveTXT(&sys); */

    /* for (int i = 0; i < sys.N_pbc; i++){
        if (sys.sites[i].index != 0) continue;
        jcv_graphedge* graph_edge = sys.sites[i].edges;
        printf("np.array([");
        while (graph_edge){
            printf("[%.16lf, %.16lf],\n", graph_edge->pos[0].x, graph_edge->pos[0].y);
            graph_edge = graph_edge->next;
        }
        printf("])");
        printf(" area = %.16lf", jcv_area(sys.sites + i));
        printf(" perimeter = %.16lf\n", jcv_perimeter(sys.sites + i));
    } */
	fclose(sys.info_snapshot.file);
	jcv_diagram_free(sys.diagram);
	return 0;
}


void constantInit(int argc, char *argv[], data* sys){


    struct option longopt[] = {
        {"number", required_argument, NULL, 'N'},
        {"Po", required_argument, NULL, 'p'},
        {"Ao", required_argument, NULL, 'a'},
        {"Ka", required_argument, NULL, 'K'},
        {"Kp", required_argument, NULL, 'P'},
        {"dt", required_argument, NULL, 'D'},
        {"time", required_argument, NULL, 'M'},
        {"dL", required_argument, NULL, 'd'},
        {"temperature", required_argument, NULL, 'T'},
        {"gamma_max", required_argument, NULL, 'g'},
        {"gamma_rate", required_argument, NULL, 'G'},
        {"version", required_argument, NULL, 'v'},
        {"k", required_argument, NULL, 'k'},
        {NULL, 0, NULL, 0}
    };
    
    int version = 1;
    int c;
    int control_dL = 0;
    while ((c = getopt_long(argc, argv, "N:q:K:P:D:m:d:T:g:G:v:k:", longopt, NULL)) != -1) {
        switch (c) {
            case 'N':
                sscanf(optarg, "%d", &sys->N);
                break;
            case 'q':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.qo);
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &sys->parameter.qo);
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &sys->parameter.qo);
                #elif JCV_type == 3
                    sys->parameter.qo = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'K':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.Ka);
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &sys->parameter.Ka);
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &sys->parameter.Ka);
                #elif JCV_type == 3
                    sys->parameter.Ka = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'P':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.Kp);
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &sys->parameter.Kp);
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &sys->parameter.Kp);
                #elif JCV_type == 3
                    sys->parameter.Kp = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'D':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->dt);
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &sys->dt);
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &sys->dt);
                #elif JCV_type == 3
                    sys->dt = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'm':
                sscanf(optarg, "%d", &sys->M);
                break;
            case 'd':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->dL);
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &sys->dL);
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &sys->dL);
                #elif JCV_type == 3
                    sys->dL = strtoflt128(optarg, NULL);
                #endif
                control_dL = 1;
                break;
            case 'T':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->parameter.T));
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &(sys->parameter.T));
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &(sys->parameter.T));
                #elif JCV_type == 3
                    sys->parameter.T = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'g':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->gamma_max));
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &(sys->gamma_max));
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &(sys->gamma_max));
                #elif JCV_type == 3
                    sys->gamma_max = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'G':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->parameter.gamma_rate));
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &(sys->parameter.gamma_rate));
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &(sys->parameter.gamma_rate));
                #elif JCV_type == 3
                    sys->parameter.gamma_rate = strtoflt128(optarg, NULL);
                #endif
                break;
            case 'v':
                sscanf(optarg, "%d", &version);
                init_genrand(version);
                break;
            case 'k':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->parameter.k));
                #elif JCV_type == 1
                    sscanf(optarg, "%lf", &(sys->parameter.k));
                #elif JCV_type == 2
                    sscanf(optarg, "%Lf", &(sys->parameter.k));
                #elif JCV_type == 3
                    sys->parameter.k = strtoflt128(optarg, NULL);
                #endif
                break;
            default:
                printf("Usage: %s [options]\n", argv[0]);
        }
    }
    sys->L = JCV_SQRT((JCV_REAL_TYPE)sys->N);
    sys->i = 0;
    sys->gamma = 0;

	int info = system("mkdir -p dump");
	if (info != 0) {
        fprintf(stderr, "Error: Failed to create directory 'dump'\n");
        return;
    }
    if (control_dL == 0){
        sys->dL = fminf(4/sys->L, 1.0);
    }

    if (sys->parameter.gamma_rate == 0.0){
        sys->shear_start = INT32_MAX;
    }
    else{
        sys->shear_start /= sys->dt;
        sys->n_to_shear_max = round(sys->gamma_max/(sys->dt*sys->parameter.gamma_rate));
        sys->parameter.gamma_rate = sys->gamma_max/(sys->dt*sys->n_to_shear_max); //so that it's exactly gamma_max at the end
        if (sys->shear_cycle == 1)
            sys->parameter.gamma_rate *= -1; //because I'm reversing the shear at shear_start in integrator.c
    }
    if (sys->info_snapshot.n_start == -1){
        sys->info_snapshot.n_start = sys->shear_start;
    }

    if (sys->info_snapshot.n_log < 0){
        sys->info_snapshot.n_log = -4*sys->n_to_shear_max/sys->info_snapshot.n_log; //need 2*n_to_shear_max to get the full cycle
    }

    if (sys->info_thermo.n_start == -1){
        sys->info_thermo.n_start = sys->shear_start;
    }

    if (sys->info_thermo.n_log == -1){
        sys->info_thermo.n_log = 4*sys->n_to_shear_max;
    }

    if (sys->shear_cycle == 0 && sys->shear_start != INT32_MAX){
        sys->M = sys->shear_start + sys->n_to_shear_max;
    }
    else if (sys->shear_cycle == 1 && sys->shear_start != INT32_MAX){
        sys->info_strobo.n_log = 4*sys->n_to_shear_max;
        sys->info_strobo.n_start = sys->shear_start;
    }
    else if (sys->shear_start == INT32_MAX){
        sys->info_strobo.n_log = INT32_MAX;
        sys->info_strobo.n_start = INT32_MAX;
    }

	char filename[200];
    do {
        sprintf(filename, "dump/N_%dqo_%Lfgamma_%Lfgammarate_%LfKa_%Lfv_%d", sys->N, (long double)sys->parameter.qo, (long double)sys->gamma_max, (long double)jcv_abs(sys->parameter.gamma_rate), (long double)sys->parameter.Ka, version);
        char snapshot_filename[256], thermo_filename[256], strobo_filename[256];
        sprintf(snapshot_filename, "%s.dump", filename);
        sprintf(thermo_filename, "%s.thermo", filename);
        sprintf(strobo_filename, "%s.strob", filename);

        sys->info_snapshot.file = fopen(snapshot_filename, "wx");
        if (sys->info_snapshot.file == NULL && errno == EEXIST) {
            version++;
        } else {
            sys->info_thermo.file = fopen(thermo_filename, "wx");
            sys->info_strobo.file = fopen(strobo_filename, "wx");
            break;
        }
    } while (1);

    fprintf(sys->info_thermo.file, "i t E gamma_actual ");
    if (sys->info_thermo.compute_stress){
        fprintf(sys->info_thermo.file, "P shear ");
    }
    if (sys->info_thermo.compute_dist_travelled){
        fprintf(sys->info_thermo.file, "frac_active ");
    }
    fprintf(sys->info_thermo.file, "\n");

    fprintf(sys->info_strobo.file, "i E gamma_actual ");
    if (sys->info_strobo.compute_stress){
        fprintf(sys->info_strobo.file, "P shear ");
    }
    if (sys->info_strobo.compute_dist_travelled){
        fprintf(sys->info_strobo.file, "frac_active ");
    }
    fprintf(sys->info_strobo.file, "\n");
   
}

