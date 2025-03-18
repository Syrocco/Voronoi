#include "force.h"
#include "helper.h"
#include "jc_voronoi.h"
#include <stdio.h>
#include "voronoi.h"
#include "thermo.h"
#include "logger.h"

jcv_point force_h(const jcv_real A, const jcv_real Aj, const jcv_real P, const jcv_real Pj, const jcv_point* h7, const jcv_point* h2, const jcv_point* h3, const parameter* param){
    
    jcv_real h72 = jcv_point_dist(h2, h7);
    jcv_real h23 = jcv_point_dist(h2, h3);
    if (jcv_real_eq(h72, 0) || jcv_real_eq(h23, 0)) { 
        printf("reset force due to vertice proximity\n");
        return (jcv_point){0, 0}; //Hacky, can be exactly 0 due to floating point stupidity
    } 
    jcv_point opp = {h3->y - h7->y, -h3->x + h7->x};
    jcv_point area = jcv_mul(param->Ka*(A - Aj), opp);

    jcv_point perimeter_last = jcv_add(jcv_mul(1/h72, jcv_sub(*h2, *h7)), jcv_mul(1/h23, jcv_sub(*h2, *h3)));

    jcv_point perimeter = jcv_mul(2*param->Kp*(P - Pj), perimeter_last);
    return jcv_add(area, perimeter);
}

jcv_point get_edge_force_ji(const jcv_site* si, const jcv_graphedge* edgei, const jcv_real Aj, const parameter* param){

    //edgei is owned by site i: si

    jcv_site* sj = edgei->neighbor;
    jcv_graphedge* edgej = sj->edges;

    jcv_graphedge* prev_edge = NULL;
    jcv_graphedge* next_edge = NULL;

    jcv_graphedge* runner = edgej;

    while (runner != NULL) {
        if (runner->neighbor == si) {
            runner = runner->next;
            break;
        }
        prev_edge = runner;
        runner = runner->next;
    }

    // If runner is not NULL, next_edge is the runner. 
    // Otherwise it's the edge from which we started
    if (runner != NULL) {
        next_edge = runner;
    }
    else{
        next_edge = sj->edges;
    }

    // If prev_edge is still NULL, it means the edge
    // in common between i and j was the first one, 
    // and prev_edge is the last one.
    if (prev_edge == NULL) {
        while (runner->next != NULL) {
            runner = runner->next;
        }
        prev_edge = runner;
    }

    const jcv_point* h7 = &(prev_edge->pos[0]);
    const jcv_point* h2 = &(prev_edge->pos[1]);
    const jcv_point* h3 = &(next_edge->pos[0]);
    const jcv_point* h8 = &(next_edge->pos[1]);

    const jcv_point* ri = &(si->p);
    const jcv_point* rj = &(sj->p);
    const jcv_point* rl = &(prev_edge->neighbor->p);
    const jcv_point* rk = &(next_edge->neighbor->p);

    const jcv_real A = jcv_area(sj);
    const jcv_real P = jcv_perimeter(sj);

    jcv_point dEj_dh2 = force_h(A, Aj, P, param->qo*JCV_SQRT(Aj), h7, h2, h3, param);
    jcv_point dEj_dh3 = force_h(A, Aj, P, param->qo*JCV_SQRT(Aj), h2, h3, h8, param);
    
    jcv_point jacobian_h2[2];
    jcv_point jacobian_h3[2];

    derivative(ri, rl, rj, jacobian_h2);
    derivative(ri, rj, rk, jacobian_h3);

    jcv_real dEj_drix = jcv_dot(dEj_dh2, jacobian_h2[0]) + jcv_dot(dEj_dh3, jacobian_h3[0]);
    jcv_real dEj_driy = jcv_dot(dEj_dh2, jacobian_h2[1]) + jcv_dot(dEj_dh3, jacobian_h3[1]);
    

    return (jcv_point){-dEj_drix, -dEj_driy};

}

jcv_point get_edge_force_ii(const jcv_site* si, const jcv_real Ai, const parameter* param){
    
    const jcv_point* ri = &(si->p);
    jcv_point* rj = NULL;
    jcv_point* rk = NULL;

    jcv_graphedge* edge = si->edges;
    jcv_graphedge* edge_after = edge->next;


    const jcv_real A = jcv_area(si);
    const jcv_real P = jcv_perimeter(si);

    jcv_real dE_drx = 0;
    jcv_real dE_dry = 0;
    jcv_point jacobian[2];

    while (edge_after != NULL){

        rj = &(edge->neighbor->p);
        rk = &(edge_after->neighbor->p);
        //printf("g = (%.16lf, %.16lf)\n", edge->pos[1].x, edge->pos[1].y);
        derivative(ri, rj, rk, jacobian);

        jcv_point dE_dh = force_h(A, Ai, P, param->qo*JCV_SQRT(Ai), &(edge->pos[0]), &(edge->pos[1]), &(edge_after->pos[1]), param);

        dE_drx += jcv_dot(dE_dh, jacobian[0]);
        dE_dry += jcv_dot(dE_dh, jacobian[1]);

        edge = edge_after;
        edge_after = edge_after->next;
    }

    //Last triangle in Delaunay made of last edge and first edge
    derivative(ri, rk, &(si->edges->neighbor->p), jacobian);
    jcv_point dE_dh = force_h(A, Ai, P, param->qo*JCV_SQRT(Ai), &(edge->pos[0]), &(edge->pos[1]), &(si->edges->pos[1]), param);
    dE_drx += jcv_dot(dE_dh, jacobian[0]);
    dE_dry += jcv_dot(dE_dh, jacobian[1]);
   
    return (jcv_point){-dE_drx, -dE_dry};

}
void derivative(const jcv_point* ri, const jcv_point* rj, const jcv_point* rk, jcv_point jacobian[2]){
    jcv_point rij = jcv_sub(*ri, *rj);
    jcv_point rik = jcv_sub(*ri, *rk);
    jcv_point rjk = jcv_sub(*rj, *rk);
    jcv_point rji = jcv_sub(*rj, *ri);
    jcv_point rkj = jcv_sub(*rk, *rj);
    jcv_point rki = jcv_sub(*rk, *ri);
    jcv_real rij_cross_rjk = jcv_cross(rij, rjk);
    jcv_real rjk_sq = jcv_lenght_sq(&rjk);
    jcv_real rik_sq = jcv_lenght_sq(&rik);
    jcv_real rij_sq = jcv_lenght_sq(&rij);
    jcv_real D = 2*rij_cross_rjk*rij_cross_rjk;

    if (D == 0){
        printf("\nri = (%f, %f), rj =  (%f, %f), rk = (%f, %f)\n", ri->x, ri->y, rj->x, rj->y, rk->x, rk->y);
        exit(3);
    }

    jcv_real alpha = rjk_sq*jcv_dot(rij, rik)/D;
    jcv_real beta = rik_sq*jcv_dot(rji, rjk)/D;
    jcv_real gamma = rij_sq*jcv_dot(rki, rkj)/D;

    jcv_point dDdri = {4*rjk.y*rij_cross_rjk, -4*rjk.x*rij_cross_rjk}; 

    jcv_point dalphadri = jcv_sub(jcv_mul(rjk_sq/D, jcv_add(rij, rik)), jcv_mul(alpha/D, dDdri));
    jcv_point dbetadri = jcv_sub(jcv_add(jcv_mul(2*jcv_dot(rji, rjk)/D, rik), jcv_mul(-rik_sq/D, rjk)), jcv_mul(beta/D, dDdri));
    jcv_point dgammadri = jcv_sub(jcv_add(jcv_mul(2*jcv_dot(rki, rkj)/D, rij), jcv_mul(-rij_sq/D, rkj)), jcv_mul(gamma/D, dDdri));

    

    // h = alpha*ri + beta*rj + gamma*rk
    // first index = derivative, .x or .y => h.x or h.y
    // ∂h.x/∂ri.x
    jacobian[0].x = alpha + dalphadri.x*ri->x + dbetadri.x*rj->x + dgammadri.x*rk->x;
    // ∂h.y/∂ri.x  
    jacobian[0].y = dalphadri.x*ri->y + dbetadri.x*rj->y + dgammadri.x*rk->y;
    // ∂h.x/∂ri.y 
    jacobian[1].x = dalphadri.y*ri->x + dbetadri.y*rj->x + dgammadri.y*rk->x;
    // ∂h.y/∂ri.y 
    jacobian[1].y = alpha + dalphadri.y*ri->y + dbetadri.y*rj->y + dgammadri.y*rk->y;
    //printf("h = (%.16lf, %.16lf)\n", alpha*ri->x + beta*rj->x + gamma*rk->x, alpha*ri->y + beta*rj->y + gamma*rk->y);

}

void compute_force(data* sys){
    #pragma omp parallel for schedule(dynamic) proc_bind(close)
    for (int i = 0; i < sys->N_pbc; i++){
        int index = sys->sites[i].index;
        if (index >= sys->N) continue;
        sys->forces[index] = get_edge_force_ii(sys->sites + i, sys->prefered_area[index], &sys->parameter); //∂Ei/∂ri
        jcv_graphedge* graph_edge = sys->sites[i].edges;
        while (graph_edge){ //looping over neighbors: graph_edge is in common between i and j
            jcv_point force_ji = get_edge_force_ji(sys->sites + i, graph_edge, sys->prefered_area[graph_edge->neighbor->index], &sys->parameter); //∂Ej/∂ri
            sys->forces[index].x += force_ji.x;                                   
            sys->forces[index].y += force_ji.y;
            graph_edge = graph_edge->next;
        }
    }
}

void check_force(data* sys){
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, 0, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    compute_force(sys);
    jcv_real old_energy = energy_total(sys)/sys->N;
    jcv_real dx = 1e-7;
    jcv_real dy = 1e-7;

    for (int i = 0; i < sys->N; i++){
        sys->positions[i].x += dx;
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++){
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, 0, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        jcv_real new_energy = energy_total(sys)/sys->N; 
        sys->positions[i].x -= dx;
        jcv_real dE_dx = sys->N*(new_energy - old_energy)/dx;

        sys->positions[i].y += dy;
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++){
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, 0, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        new_energy = energy_total(sys)/sys->N; 
        sys->positions[i].y -= dy;
        jcv_real dE_dy =  sys->N*(new_energy - old_energy)/dy;

        if (jcv_abs((dE_dx + sys->forces[i].x)/dE_dx) > 1e-2 && jcv_abs((dE_dy +sys->forces[i].y)/dE_dy) > 2e-1){
            printf("Force is wrong for %d\n", i);
            printf("dE_dx = %.14lf, force_x = %.14lf\n", -dE_dx, sys->forces[i].x);
            printf("dE_dy = %.14lf, force_y = %.14lf\n", -dE_dy, sys->forces[i].y);
        }
    }
    saveTXT(sys);
}

