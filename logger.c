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

    if ((sys->i - sys->info_strobo.n_start)%(sys->info_strobo.n_log) == 0 
         && sys->i >= sys->info_strobo.n_start){
        computeStrobo(sys);
    }

}

void computeThermo(data* sys){
    
    jcv_real E = energy_total(sys)/sys->N;
    fprintf(sys->info_thermo.file, "%d %Lf %.16Lf %Lf ", sys->i, (long double)sys->time, (long double)E, (long double)sys->gamma);
    printf("i = %d, dt = %Lf, E = %.16Lf, γ = %.16Lf ", sys->i, (long double)sys->dt, (long double)E, (long double)sys->gamma);

    if (sys->info_thermo.compute_stress){
        jcv_real stress[2][2];
        stress_total(sys, stress);
        jcv_real shear_stress = (stress[0][1] + stress[1][0])/2;
        jcv_real pressure = (stress[0][0] + stress[1][1])/2;
        fprintf(sys->info_thermo.file, "%.16Lf %.16Lf ", (long double)pressure, (long double)shear_stress);
        printf("P = %Lf, shear = %Lf ", (long double)pressure, (long double)shear_stress);
    }
    
    if (sys->info_thermo.compute_dist_travelled){
        jcv_real frac_active = 0;
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
        frac_active = frac_active/sys->N;
        fprintf(sys->info_thermo.file, "%Lf ", (long double)frac_active);
        printf("frac_active = %Lf ", (long double)frac_active);
    }
    fprintf(sys->info_thermo.file, "\n");
    printf("\n");
    fflush(sys->info_thermo.file);    
}

void computeStrobo(data* sys){
    
    jcv_real E = energy_total(sys)/sys->N;
    fprintf(sys->info_strobo.file, "%d %Lf %Lf ", sys->i, (long double)E, (long double)sys->gamma);

    if (sys->info_strobo.compute_stress){
        jcv_real stress[2][2];
        stress_total(sys, stress);
        jcv_real shear_stress = (stress[0][1] + stress[1][0])/2;
        jcv_real pressure = (stress[0][0] + stress[1][1])/2;
        fprintf(sys->info_strobo.file, "%Lf %Lf ", (long double)pressure, (long double)shear_stress);
    }
    
    if (sys->info_strobo.compute_dist_travelled){
        jcv_real frac_active = 0;
        jcv_real distance_move[sys->N];
        distance_moved(sys, sys->old_info.old_positions_strobo, distance_move);
        for (int i = 0; i < sys->N; i++){
            if (distance_move[i] > 0.001){
                frac_active++;
            }
        }
        frac_active = frac_active/sys->N;
        fprintf(sys->info_strobo.file, "%Lf ", (long double)frac_active);
        if (frac_active == 0){
            printf("No more active particles");
            exit(3);
        } 
    }
    fprintf(sys->info_strobo.file, "\n");
    fflush(sys->info_strobo.file);    
}

void saveTXT(data* sys){
    FILE* file = sys->info_snapshot.file;
    jcv_real gamma = sys->gamma;
    jcv_real L = sys->L;    
    jcv_point* p = sys->positions;
    jcv_point* f = sys->forces;
    jcv_real* a = sys->prefered_area;
    int N = sys->N;
    int m = sys->i;
    int N_pbc = sys->N_pbc;
    int NN = {sys->info_snapshot.include_boundary ? N_pbc : N};
    jcv_real dL = gamma*L;

    jcv_real stress[N][2][2];
    jcv_real distance_move[N];
    jcv_real area[N];
    jcv_real perimeter[N];

 

    fprintf(file, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz\n%.16Lf %.16Lf %.16Lf\n0 %.16Lf 0\n0 0 0\nITEM: ATOMS id prefered_area x y fx fy ", m, NN, (long double)jcv_min(0, dL), (long double)jcv_max(L, L + dL), (long double)dL, (long double)L);
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
            fprintf(file, "%d %Lf %Lf %Lf nan nan nan nan ", i, (long double)a[i], (long double)p[i].x, (long double)p[i].y);
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
            fprintf(file, "%d %.16Lf %.16Lf %.16Lf %.16Lf %.16Lf ", i, (long double)a[i], (long double)p[i].x, (long double)p[i].y, (long double)f[i].x, (long double)f[i].y);
            if (sys->info_snapshot.compute_stress){
                fprintf(file, "%.16Lf ", (long double)(stress[i][0][1] + stress[i][1][0])/2);
            }
            if (sys->info_snapshot.compute_dist_travelled){
                fprintf(file, "%Lf ", (long double)distance_move[i]);
            }
            if (sys->info_snapshot.compute_area){
                fprintf(file, "%.16Lf ", (long double)area[i]);
            }
            if (sys->info_snapshot.compute_perimeter){
                fprintf(file, "%.16Lf ", (long double)perimeter[i]);
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
        fprintf(file, "%Lf %Lf\n", (long double)points[i].x, (long double)points[i].y);
    }

    const jcv_graphedge* graph_edge;
    for (int i = 0; i < N_pbc; i++){
        //if (sites[i].index >= N) continue;

        graph_edge = sites[i].edges;
        while (graph_edge) {
            fprintf(file, "%Lf %Lf\n", (long double)graph_edge->pos[0].x, (long double)graph_edge->pos[0].y);
            fprintf(file, "%Lf %Lf\n", (long double)graph_edge->pos[1].x, (long double)graph_edge->pos[1].y);
            graph_edge = graph_edge->next;
        }
    }
    fclose(file);
}

void unexpectedClosure(data* sys){
    fprintf(sys->info_thermo.file, "nan nan nan nan ");
    if (sys->info_thermo.compute_stress){
        fprintf(sys->info_thermo.file, "nan nan ");
    }
    
    if (sys->info_thermo.compute_dist_travelled){
        fprintf(sys->info_thermo.file, "nan ");
    }
    saveTXT(sys);
}