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



int main(int argc, char *argv[]){
    init_genrand(0);

    data sys;
    sys.parameter.Ao = 1.0;
    sys.parameter.Po = 4;
    sys.parameter.Ka = 1;
    sys.parameter.Kp = 1;
    sys.N = 100;
    sys.M = 620; 
    sys.dt = 0.005;
    sys.T = 0.;
    sys.gamma_rate = 0.0;
    
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
    


    //randomInitial(&sys);
    
    rsaInitial(&sys, 0.4);
    

    for (sys.i = 0; sys.i < sys.M; sys.i++){
        
        
        
        eulerStep(&sys);
        //fireStep(&sys);

        
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
        {"temperature", required_argument, NULL, 'T'},
        {NULL, 0, NULL, 0}
    };
    

    int c;
    int control_dL = 0;
    while ((c = getopt_long(argc, argv, "N:p:a:k:P:D:m:d:T:", longopt, NULL)) != -1){
        switch(c){
            case 'N':
                sscanf(optarg, "%d", &sys->N);
                break;
            case 'p':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.Po);
                #else
                    sscanf(optarg, "%lf", &sys->parameter.Po);
                #endif
                break;
            case 'a':
                #if JCV_type == 0
                    sscanf(optarg, "%f", &sys->parameter.Ao);    
                #else
                    sscanf(optarg, "%lf", &sys->parameter.Ao);
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
                    sscanf(optarg, "%f", &sys->T);
                #else
                    sscanf(optarg, "%lf", &sys->T);
                #endif
                break;
        }
    }
    sys->L = JCV_SQRT((JCV_REAL_TYPE)sys->N);
    sys->i = 0;
    sys->amount_of_def = 0;

	int info = system("mkdir -p dump");
	if (info != 0) {
        fprintf(stderr, "Error: Failed to create directory 'dump'\n");
        return;
    }
    if (control_dL == 0){
        sys->dL = fminf(5/sys->L, 1.0);
    }
	char filename[100];
	sprintf(filename, "dump/N_%dPo_%f.dump", sys->N, sys->parameter.Po);
	sys->file = fopen(filename, "w");
}

