#include "voronoi.h"
#include "logger.h"
#include "force.h"

void loggers(data* sys){
    //only call saveTXT every sn_log steps beginning at n_start
    if ((sys->i - sys->info_particles.n_start)%(sys->info_particles.n_log) == 0 
         && sys->i >= sys->info_particles.n_start){
        saveTXT(sys);
    }

    if ((sys->i - sys->info_thermo.n_start)%(sys->info_thermo.n_log) == 0 
         && sys->i >= sys->info_thermo.n_start){
        computeThermo(sys);
    }
}

void computeThermo(data* sys){
    jcv_real E = energy_total(sys->sites, sys->N, sys->N_pbc, &sys->parameter);
    jcv_real stress[2][2];
    stress_total(sys, stress);
    jcv_real shear_stress = (stress[0][1] + stress[1][0])/2;
    jcv_real pressure = (stress[0][0] + stress[1][1])/2;
    fprintf(sys->info_thermo.file, "%d %lf %lf %lf %lf\n", sys->i, E, pressure, shear_stress, sys->gamma);
    fflush(sys->info_thermo.file);

    printf("i = %d, E = %f, P = %f, shear = %f, gamma = %f\n", sys->i, E, pressure, shear_stress, sys->gamma);
}

void saveTXT(data* sys){
    FILE* file = sys->info_particles.file;
    jcv_real gamma = sys->gamma;
    jcv_real L = sys->L;    
    jcv_point* p = sys->positions;
    jcv_point* f = sys->forces;
    jcv_point* v = sys->velocities;
    int N = sys->N;
    int m = sys->i;
    int N_pbc = sys->N_pbc;
    int NN = {sys->info_particles.include_boundary ? N_pbc : N};
    jcv_real dL = gamma*L;

    jcv_real stress[N][2][2];

    for (int i = 0; i < N_pbc; i++){
        if (sys->sites[i].index >= N) continue;
        stress_unique(sys->sites + i, &sys->parameter, stress[sys->sites[i].index]);
    }
    fprintf(file, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz\n0 %f %f\n0 %f 0\n0 0 0\nITEM: ATOMS id x y vx vy fx fy shear\n", m, NN, L + dL, dL, L);
    for(int i = 0; i < NN; i++){
        if (i >= N){
            fprintf(file, "%d %lf %lf nan nan nan nan nan\n", i, p[i].x, p[i].y);
            continue;
        }
        else{
            fprintf(file, "%d %lf %lf %lf %lf %lf %lf %lf\n", i, p[i].x, p[i].y, v[i].x, v[i].y, f[i].x, f[i].y, (stress[i][0][1] + stress[i][1][0])/2);
        }
    }
    fflush(file);
}

void write(FILE* file, const char* filename, const jcv_point* points, const jcv_site* sites, int N, int N_pbc){
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }

    fprintf(file, "%d %d\n", N, N_pbc);
    for (int i = 0; i < N_pbc; i++) {
        fprintf(file, "%f %f\n", (float)points[i].x, (float)points[i].y);
    }

    const jcv_graphedge* graph_edge;
    for (int i = 0; i < N_pbc; i++){
        //if (sites[i].index >= N) continue;

        graph_edge = sites[i].edges;
        while (graph_edge) {
            fprintf(file, "%f %f\n", (float)graph_edge->pos[0].x, (float)graph_edge->pos[0].y);
            fprintf(file, "%f %f\n", (float)graph_edge->pos[1].x, (float)graph_edge->pos[1].y);
            graph_edge = graph_edge->next;
        }
    }
    fclose(file);
}