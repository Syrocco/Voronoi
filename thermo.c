#include "thermo.h"
#include "voronoi.h"
#include "jc_voronoi.h"
#include "helper.h"

void shear(data* sys){
    static int start = 1;
    static int lag = 0;
    int normalized_time = sys->i - sys->shear_start - lag;
    if (normalized_time >= 0){
        if (normalized_time%sys->n_to_shear_max == 0){
            sys->parameter.gamma_rate *= -1;
        }
        sys->gamma += sys->dt*sys->parameter.gamma_rate;
    }
    if (start && (normalized_time == sys->n_to_shear_max)){ //when we went from 0 to gamma_max the first time
        lag = sys->n_to_shear_max;
        sys->n_to_shear_max *= 2;
        start = 0;
    }
}

inline jcv_real energy_unique(const jcv_site* site, const jcv_real A, const parameter* param){
    jcv_real area = jcv_area(site);
    jcv_real perimeter = jcv_perimeter(site);
    return param->Ka*(area - A)*(area - A)
           + param->Kp*(perimeter - param->qo*JCV_SQRT(A))*(perimeter - param->qo*JCV_SQRT(A));
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
    const jcv_site* site = sys->sites + num;
    

    jcv_real Ai = jcv_area(site);
    jcv_real Pi = jcv_perimeter(site);
    jcv_real pressure = sys->parameter.Ka*(Ai - sys->prefered_area[site->index]);
    stress[0][0] = pressure;
    stress[1][1] = pressure;

    stress[0][1] = 0;
    stress[1][0] = 0;
    jcv_real prefac = 1/(2*Ai); // 2 because each edge is counted twice
    jcv_graphedge* edge = site->edges;
    while (edge){

        jcv_point l_ij = jcv_sub(edge->pos[1], edge->pos[0]);
        jcv_point l_ij_normalized = jcv_mul(1/jcv_lenght(&l_ij), l_ij);
        jcv_real Pj = jcv_perimeter(edge->neighbor);
        jcv_point tension_ij = jcv_mul(2*sys->parameter.Kp*((Pj - sys->parameter.qo*JCV_SQRT(sys->prefered_area[edge->neighbor->index]))
                                                          + (Pi - sys->parameter.qo*JCV_SQRT(sys->prefered_area[site->index]))), l_ij_normalized);
       
        stress[0][0] += prefac*tension_ij.x*l_ij.x;
        stress[0][1] += prefac*tension_ij.x*l_ij.y;
        stress[1][0] += prefac*tension_ij.y*l_ij.x;    
        stress[1][1] += prefac*tension_ij.y*l_ij.y;

        edge = edge->next;
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
        jcv_real area = jcv_area(sys->sites + i);
        stress_unique(sys, i, temp_stress);
        stress[0][0] += area*temp_stress[0][0];
        stress[0][1] += area*temp_stress[0][1];
        stress[1][0] += area*temp_stress[1][0];
        stress[1][1] += area*temp_stress[1][1];
    }
    stress[0][0] /= sys->N; //N is total area in reduced units.
    stress[0][1] /= sys->N;
    stress[1][0] /= sys->N;
    stress[1][1] /= sys->N;
}

void distance_moved(data* sys, jcv_point* old_positions, jcv_real* dist){
    for (int i = 0; i < sys->N; i++){
        jcv_point diff = jcv_sub(sys->positions[i], old_positions[i]);
        pbc_distance(&diff, sys->L);
        dist[i] = jcv_lenght(&diff);
        
        old_positions[i] = sys->positions[i];
    }
}