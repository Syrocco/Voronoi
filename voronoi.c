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
jcv_real energyUnique(data* sys, const parameter* param, int num_particle);
int main(int argc, char *argv[]){
    init_genrand(0);

    data sys;
    sys.parameter.qo = 3.65;
    sys.parameter.Ka = 1;
    sys.parameter.Kp = 1;

    sys.size_large_over_small = 4.0/3.0;
    sys.n_frac_small = 0.5;

    sys.N = 300;
    sys.M = 100000; 
    sys.dt = 0.001;
    sys.parameter.T = 0.0;
    sys.parameter.gamma_rate = 0.0;
    sys.gamma_max = 0.2;
    sys.shear_start = 0;
    sys.dt_fire = 0.1;
    sys.shear_cycle = 1;

    sys.info_snapshot.n_log = -1;
    sys.info_snapshot.n_start = -1;
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
    rsaInitial(&sys, 0.4);
    //read_from_dump_initial(&sys, "N_300qo_3.550000gamma_3.000000gammarate_0.100000Ka_3.000000v_2.dump", 35);

    for (sys.i = 0; sys.i < sys.M; sys.i++){
        
        //backwardEulerStep(&sys);
        //rkf45Step(&sys);
        //rk4Step(&sys);
        //eulerStep(&sys);
        //fireStep(&sys);
        conjugateGradientStep(&sys);
        
    }
    
    check_force(&sys);

	fclose(sys.info_snapshot.file);
	jcv_diagram_free(sys.diagram);
	return 0;
}


void constantInit(int argc, char *argv[], data* sys){


    struct option longopt[] = {
        {"number", required_argument, NULL, 'N'},
        {"Po", required_argument, NULL, 'p'},
        {"Ao", required_argument, NULL, 'a'},
        {"Ka", required_argument, NULL, 'k'},
        {"Kp", required_argument, NULL, 'P'},
        {"dt", required_argument, NULL, 'D'},
        {"time", required_argument, NULL, 'M'},
        {"dL", required_argument, NULL, 'd'},
        {"temperature", required_argument, NULL, 'T'},
        {"gamma_max", required_argument, NULL, 'g'},
        {NULL, 0, NULL, 0}
    };
    
    int version = 1;
    int c;
    int control_dL = 0;
    while ((c = getopt_long(argc, argv, "N:q:k:P:D:m:d:T:g:G:v:", longopt, NULL)) != -1){
        switch(c){
            case 'N':
                sscanf(optarg, "%d", &sys->N);
                break;
            case 'q':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.qo);
                #else
                    sscanf(optarg, "%lf", &sys->parameter.qo);
                #endif
                break;
            case 'k':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.Ka);
                #else
                    sscanf(optarg, "%lf", &sys->parameter.Ka);
                #endif
                break;
            case 'P':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.Kp);
                #else
                    sscanf(optarg, "%lf", &sys->parameter.Kp);
                #endif
                break;
            case 'D':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->dt);
                #else
                    sscanf(optarg, "%lf", &sys->dt);
                #endif
                break;
            case 'm':
                sscanf(optarg, "%d", &sys->M);
                break;
            case 'd':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->dL);
                #else
                    sscanf(optarg, "%lf", &sys->dL);
                #endif
                control_dL = 1;
                break;
            case 'T':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->parameter.T));
                #else
                    sscanf(optarg, "%lf", &(sys->parameter.T));
                #endif
                break;
            case 'g':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->gamma_max));
                #else
                    sscanf(optarg, "%lf", &(sys->gamma_max));
                #endif
                break;
            case 'G':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &(sys->parameter.gamma_rate));
                #else
                    sscanf(optarg, "%lf", &(sys->parameter.gamma_rate));
                #endif
                break;
            case 'v':
                sscanf(optarg, "%d", &version);
                init_genrand(version);
                break;
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
        sprintf(filename, "dump/N_%dqo_%fgamma_%fgammarate_%lfKa_%lfv_%d", sys->N, sys->parameter.qo, sys->gamma_max, -sys->parameter.gamma_rate, sys->parameter.Ka, version);
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

