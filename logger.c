#include "voronoi.h"
#include "logger.h"
#include "force.h"
#include "thermo.h"
#include "helper.h"

void loggers(data* sys){
    //only call saveTXT every sn_log steps beginning at n_start
    if ((sys->i - sys->info_snapshot.n_start)%(sys->info_snapshot.n_log) == 0 
         && sys->i >= sys->info_snapshot.n_start){
        saveTXT(sys);
    }

    if ((sys->i - sys->info_thermo.n_start)%(sys->info_thermo.n_log) == 0 
         && sys->i >= sys->info_thermo.n_start){
        computeThermo(sys);
    }

}

void computeThermo(data* sys){
    
    

    jcv_real E = energy_total(sys);
    fprintf(sys->info_thermo.file, "%d %lf %lf ", sys->i, E, sys->gamma);
    printf("i = %d, E = %lf, γ = %lf ", sys->i, E, sys->gamma);

    if (sys->info_thermo.compute_stress){
        jcv_real stress[2][2];
        stress_total(sys, stress);
        jcv_real shear_stress = (stress[0][1] + stress[1][0])/2;
        jcv_real pressure = (stress[0][0] + stress[1][1])/2;
        fprintf(sys->info_thermo.file, "%lf %lf ", pressure, shear_stress);
        printf("P = %lf, shear = %lf ", pressure, shear_stress);
    }
    
    if (sys->info_thermo.compute_dist_travelled){
        int frac_active = 0;
        jcv_real distance_move[sys->N];
        distance_moved(sys, sys->old_info.old_positions_thermo, distance_move);
        for (int i = 0; i < sys->N; i++){
            #if JCV_type == 0
                if (distance_move[i] > 0.001){
            #else
                if (distance_move[i] > 0.0000001){
            #endif
                frac_active++;
            }
        }
        fprintf(sys->info_thermo.file, "%d ", frac_active);
        printf("frac_active = %d ", frac_active);
    }
    fprintf(sys->info_thermo.file, "\n");
    printf("\n");
    fflush(sys->info_thermo.file);

    
}

void saveTXT(data* sys){
    FILE* file = sys->info_snapshot.file;
    jcv_real gamma = sys->gamma;
    jcv_real L = sys->L;    
    jcv_point* p = sys->positions;
    jcv_point* f = sys->forces;
    jcv_point* v = sys->velocities;
    int N = sys->N;
    int m = sys->i;
    int N_pbc = sys->N_pbc;
    int NN = {sys->info_snapshot.include_boundary ? N_pbc : N};
    jcv_real dL = gamma*L;

    jcv_real stress[N][2][2];
    jcv_real distance_move[N];
    jcv_real area[N];
    jcv_real perimeter[N];

 

    fprintf(file, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz\n0 %f %f\n0 %f 0\n0 0 0\nITEM: ATOMS id x y vx vy fx fy ", m, NN, L + dL, dL, L);
    if (sys->info_snapshot.compute_stress){
        fprintf(file, "shear ");
        for (int i = 0; i < N_pbc; i++){
            if (sys->sites[i].index >= N) continue;
            stress_unique(sys, i, stress[sys->sites[i].index]);
        }
    }
    if (sys->info_snapshot.compute_dist_travelled){
        fprintf(file, "dist ");
        distance_moved(sys, sys->old_info.old_positions_snapshot, distance_move);
    }
    if (sys->info_snapshot.compute_area){
        fprintf(file, "area ");
        for (int i = 0; i < N_pbc; i++){
            if (sys->sites[i].index >= N) continue;
            area[sys->sites[i].index] = jcv_area(sys->sites + i);
        }
    }
    if (sys->info_snapshot.compute_perimeter){
        fprintf(file, "perimeter ");
        for (int i = 0; i < N_pbc; i++){
            if (sys->sites[i].index >= N) continue;
            perimeter[sys->sites[i].index] = jcv_perimeter(sys->sites + i);
        }
    }
    fprintf(file, "\n");
    
    for(int i = 0; i < NN; i++){
        if (i >= N){
            fprintf(file, "%d %lf %lf nan nan nan nan ", i, p[i].x, p[i].y);
            if (sys->info_snapshot.compute_stress){
                fprintf(file, "nan ");
            }
            if (sys->info_snapshot.compute_dist_travelled){
                fprintf(file, "nan ");
            }
            if (sys->info_snapshot.compute_area){
                fprintf(file, "nan ");
            }
            if (sys->info_snapshot.compute_perimeter){
                fprintf(file, "nan ");
            }
            fprintf(file, "\n");
        }
        else{
            fprintf(file, "%d %lf %lf %lf %lf %.12lf %.12lf ", i, p[i].x, p[i].y, v[i].x, v[i].y, f[i].x, f[i].y);
            if (sys->info_snapshot.compute_stress){
                fprintf(file, "%lf ", (stress[i][0][1] + stress[i][1][0])/2);
            }
            if (sys->info_snapshot.compute_dist_travelled){
                fprintf(file, "%lf ", distance_move[i]);
            }
            if (sys->info_snapshot.compute_area){
                fprintf(file, "%lf ", area[i]);
            }
            if (sys->info_snapshot.compute_perimeter){
                fprintf(file, "%lf ", perimeter[i]);
            }
            fprintf(file, "\n");
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