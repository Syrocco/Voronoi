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


int main(int argc, char *argv[]){
    init_genrand(0);

    data sys;
    sys.parameter.Ao = 1.0;
    sys.parameter.Po = 3.7;
    sys.parameter.Ka = 1.0;
    sys.parameter.Kp = 1.0;
    sys.N = 100;
    sys.M = 10000; 
    sys.L = JCV_SQRT((JCV_REAL_TYPE)sys.N);;
    sys.dt = 0.1;
    
    constantInit(argc, argv, &sys);
    
    jcv_diagram diagram;
    memset(&(diagram), 0, sizeof(jcv_diagram));
    jcv_point positions[9*sys.N];
    jcv_point velocities[sys.N];
    jcv_point forces[sys.N];
    sys.diagram = &diagram;
    sys.positions = positions;
    sys.velocities = velocities;
    sys.forces = forces;
    

    
    for (int i = 0; i < sys.N; i++) {
        sys.positions[i].x = drand(0, sys.L);
        sys.positions[i].y = drand(0, sys.L);
    }

    
    
    
    for (int m = 0; m < sys.M; m++){
        //clock_t start_time = clock();
        
        populate_points(sys.positions, sys.N, sys.L);
        jcv_diagram_generate(9*sys.N, sys.positions, NULL, 0, sys.diagram);
        sys.sites = jcv_diagram_get_sites(sys.diagram);
        
        //clock_t end_time = clock();
        //double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        //printf("Time for populate_points, jcv_diagram_generate, and jcv_diagram_get_sites: %f seconds\n", elapsed_time);

        if (m%25 == 0){
            saveTXT(sys.file, sys.positions, sys.N, m, sys.L);
            //sprintf(filename, "dump/%d.txt", m);
            //write(file, filename, points, sites, N);
        }

        //start_time = clock();
        compute_force(&sys);
        //end_time = clock();
        //elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        //printf("Time for compute_force: %f seconds\n", elapsed_time);

        for (int i = 0; i < sys.N; i++){
            sys.positions[i].x = pbc(sys.positions[i].x + sys.dt*sys.forces[i].x + JCV_SQRT(sys.dt)*drand(-0.1, 0.1), sys.L);
            sys.positions[i].y = pbc(sys.positions[i].y + sys.dt*sys.forces[i].y + JCV_SQRT(sys.dt)*drand(-0.1, 0.1), sys.L);
        }
        printf("m = %d, E = %f \n", m, energy(sys.sites, sys.N, &sys.parameter));
        
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
        {NULL, 0, NULL, 0}
    };
    

    int c;

    while ((c = getopt_long(argc, argv, "N:p:a:k:P:D:m:", longopt, NULL)) != -1){
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
        }
    }
    sys->L = JCV_SQRT((JCV_REAL_TYPE)sys->N);
	int info = system("mkdir -p dump");
	if (info != 0) {
        fprintf(stderr, "Error: Failed to create directory 'dump'\n");
        return;
    }

	char filename[100];
	sprintf(filename, "dump/N_%dPo_%f.dump", sys->N, sys->parameter.Po);
	sys->file = fopen(filename, "w");
}

