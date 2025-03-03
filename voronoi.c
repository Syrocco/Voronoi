#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define JC_VORONOI_IMPLEMENTATION
// If you wish to use doubles
//#define JCV_REAL_TYPE double
//#define JCV_FABS fabs
//#define JCV_ATAN2 atan2
//#define JCV_CEIL ceilP
//#define JCV_FLOOR floor
//#define JCV_FLT_MAX 1.7976931348623157E+308
#include "jc_voronoi.h"
#include "mersenne.c"
#include "helper.h"
#include "force.h"
#include "voronoi.h"
#include <getopt.h>

#define TIME_FUNCTION(func, ...) \
    do { \
        clock_t start_time = clock(); \
        func(__VA_ARGS__); \
        clock_t end_time = clock(); \
        double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC; \
        printf("Time for %s: %f seconds\n", #func, elapsed_time); \
    } while (0)

int main(int argc, char *argv[]){
    init_genrand(0);

    data sys;
    sys.parameter.Ao = 1.0;
    sys.parameter.Po = 3.7;
    sys.parameter.Ka = 1.0;
    sys.parameter.Kp = 1.0;
    sys.N = 10000;
    sys.M = 1; 
    sys.L = JCV_SQRT((JCV_REAL_TYPE)sys.N);;
    sys.dt = 0.1;
    sys.gamma_rate = 3;
    sys.i = 0;
    
    constantInit(argc, argv, &sys);
    
    jcv_diagram diagram;
    memset(&(diagram), 0, sizeof(jcv_diagram));
    jcv_point positions[8*sys.N];
    jcv_point velocities[sys.N];
    jcv_point forces[sys.N];
    sys.diagram = &diagram;
    sys.positions = positions;
    sys.velocities = velocities;
    sys.forces = forces;
    
    sys.N_pbc = sys.N;
    for (int i = 0; i < sys.N; i++) {
        sys.positions[i].x = drand(0, sys.L);
        sys.positions[i].y = drand(0, sys.L);
        addBoundary(&sys, i);
    }
    
    for (sys.i = 0; sys.i < sys.M; sys.i++){
        sys.deformation_by_lenght = (sys.i + 1)*sys.dt*sys.gamma_rate;
        TIME_FUNCTION(jcv_diagram_generate, sys.N_pbc, sys.positions, NULL, 0, sys.diagram);
        sys.sites = jcv_diagram_get_sites(sys.diagram);
        printf("m = %d, E = %f \n", sys.i, energy(sys.sites, sys.N, sys.N_pbc, &sys.parameter));

        if (sys.i%1 == 0){
            saveTXT(sys.file, sys.positions, sys.N_pbc, sys.i, sys.L);

            sys.N_pbc = sys.N;
            for (int i = 0; i < sys.N; i++){
            
                sys.positions[i].x += sys.positions[i].y*sys.deformation_by_lenght;
                pbc(&sys.positions[i], sys.L, sys.deformation_by_lenght);
                addBoundary(&sys, i);
            }
			char filename[100];
            sprintf(filename, "dump/a.txt");
			FILE* file = fopen(filename, "w");
            write(file, filename, sys.positions, sys.sites, sys.N_pbc, sys.N_pbc);
			fclose(file);
            exit(3);
            
        }

        TIME_FUNCTION(compute_force,&sys);
        

        sys.N_pbc = sys.N;
        for (int i = 0; i < sys.N; i++){
            
            sys.positions[i].x = sys.positions[i].x + sys.dt*sys.forces[i].x + JCV_SQRT(sys.dt)*drand(-0.1, 0.1);
            sys.positions[i].y = sys.positions[i].y + sys.dt*sys.forces[i].y + JCV_SQRT(sys.dt)*drand(-0.1, 0.1);
            pbc(&sys.positions[i], sys.L, sys.dt*sys.i*sys.L);
            addBoundary(&sys, i);
        }
        
        
    }

	fclose(sys.file);
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
        {NULL, 0, NULL, 0}
    };
    

    int c;
    int control_dL = 0;
    while ((c = getopt_long(argc, argv, "N:p:a:k:P:D:m:d:", longopt, NULL)) != -1){
        switch(c){
            case 'N':
                sscanf(optarg, "%d", &sys->N);
                break;
            case 'p':
                sscanf(optarg, "%f", &sys->parameter.Po);
                break;
            case 'a':
                sscanf(optarg, "%f", &sys->parameter.Ao);
                break;
            case 'k':
                sscanf(optarg, "%f", &sys->parameter.Ka);
                break;
            case 'P':
                sscanf(optarg, "%f", &sys->parameter.Kp);
                break;
            case 'D':
                sscanf(optarg, "%f", &sys->dt);
                break;
            case 'm':
                sscanf(optarg, "%d", &sys->M);
                break;
            case 'd':
                sscanf(optarg, "%f", &sys->dL);
                control_dL = 1;
                break;
        }
    }
    sys->L = JCV_SQRT((JCV_REAL_TYPE)sys->N);
	int info = system("mkdir -p dump");
	if (info != 0) {
        fprintf(stderr, "Error: Failed to create directory 'dump'\n");
        return;
    }
    if (control_dL == 0){
        sys->dL = fminf(3/sys->L, 1.0);
    }
	char filename[100];
	sprintf(filename, "dump/N_%dPo_%f.dump", sys->N, sys->parameter.Po);
	sys->file = fopen(filename, "w");
}

