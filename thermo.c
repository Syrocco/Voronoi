#include "thermo.h"
#include "voronoi.h"
#include "jc_voronoi.h"
#include "helper.h"

void shear(data* sys){
    int normalized_time = sys->i - sys->shear_start;
    if (normalized_time >= 0){
        if (normalized_time%sys->n_to_shear_max == 0){
            sys->parameter.gamma_rate *= -1;
        }
        sys->gamma += sys->dt*sys->parameter.gamma_rate;
    }
}

inline jcv_real energy_unique(const jcv_site* site, const jcv_real A, const parameter* param){
    return param->Ka*(jcv_area(site) - A)*(jcv_area(site) - A)
           + param->Kp*(jcv_perimeter(site) - param->qo*JCV_SQRT(A))*(jcv_perimeter(site) - param->qo*JCV_SQRT(A));
}

jcv_real energy_total(data* sys){
    jcv_real E = 0;
    for (int i = 0; i < sys->N_pbc; i++){
        if (sys->sites[i].index >= sys->N) continue;
        E += energy_unique(sys->sites + i, sys->prefered_area[sys->sites[i].index], &(sys->parameter));
    }
    return E;
}

// https://pmc.ncbi.nlm.nih.gov/articles/PMC11398538/pdf/nihpp-2409.04383v1.pdf
// https://www.pnas.org/doi/pdf/10.1073/pnas.1705921114
inline void stress_unique(data* sys, int num, jcv_real stress[2][2]){
    
    stress[0][1] = 0;
    stress[1][0] = 0;
    const jcv_site* site = sys->sites + num;
    jcv_real Ai = jcv_area(site);
    jcv_real Pi = jcv_perimeter(site);
    jcv_real pressure = sys->parameter.Ka*(Ai - sys->prefered_area[site->index]);
    jcv_graphedge* edge = site->edges;
    
    stress[0][0] = pressure;
    stress[1][1] = pressure;

    jcv_real prefac = 1/(2*Ai); // 2 because each edge is counted twice
    while (edge){
        jcv_point l_ij = jcv_sub(edge->pos[1], edge->pos[0]);
        jcv_point l_ij_normalized = jcv_mul(1/jcv_lenght(&l_ij), l_ij);
        jcv_real Pj = jcv_perimeter(edge->neighbor);
        jcv_point tension = jcv_mul(2*sys->parameter.Kp*((Pj - sys->parameter.qo*JCV_SQRT(sys->prefered_area[edge->neighbor->index])) + (Pi - sys->parameter.qo*JCV_SQRT(sys->prefered_area[site->index]))), l_ij_normalized);
        edge = edge->next;
        stress[0][0] += prefac*tension.x*l_ij.x;
        stress[0][1] += prefac*tension.x*l_ij.y;
        stress[1][0] += prefac*tension.y*l_ij.x;    
        stress[1][1] += prefac*tension.y*l_ij.y;
    }
}

void stress_total(data* sys, jcv_real stress[2][2]){
    stress[0][0] = 0;
    stress[0][1] = 0;
    stress[1][0] = 0;
    stress[1][1] = 0;
    jcv_real temp_stress[2][2];
    for (int i = 0; i < sys->N_pbc; i++){
        if (sys->sites[i].index >= sys->N) continue;
        stress_unique(sys, i, temp_stress);
        stress[0][0] += temp_stress[0][0];
        stress[0][1] += temp_stress[0][1];
        stress[1][0] += temp_stress[1][0];
        stress[1][1] += temp_stress[1][1];
    }
}

void distance_moved(data* sys, jcv_point* old_positions, jcv_real* dist){
    for (int i = 0; i < sys->N; i++){
        jcv_point diff = jcv_sub(sys->positions[i], old_positions[i]);
        pbc_distance(&diff, sys->L);
        dist[i] = jcv_lenght(&diff);
        
        old_positions[i] = sys->positions[i];
    }
}