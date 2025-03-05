#include "force.h"
#include "helper.h"
#include "jc_voronoi.h"
#include <stdio.h>
#include "voronoi.h"

jcv_point force_h(const jcv_real A, const jcv_real P, const jcv_point* h7, const jcv_point* h2, const jcv_point* h3, const parameter* param){
    jcv_real h72 = jcv_point_dist(h2, h7);
    jcv_real h23 = jcv_point_dist(h2, h3);
    if (jcv_real_eq(h72, 0) || jcv_real_eq(h23, 0)) { 
        //return (jcv_point){0, 0}; //Hacky, can be exactly 0 due to floating point stupidity
    }
    jcv_point opp = {h3->y - h7->y, h3->x - h7->x};
    jcv_point area = jcv_mul(param->Ka*(A - param->Ao), opp);

    jcv_point perimeter_last = jcv_add(jcv_mul(1/h72, jcv_sub(*h2, *h7)), jcv_mul(1/h23, jcv_sub(*h2, *h3)));

    jcv_point perimeter = jcv_mul(2*param->Kp*(P - param->Po), perimeter_last);
    return jcv_add(area, perimeter);
}

jcv_point get_edge_force_ji(const jcv_site* si, const jcv_graphedge* edgei, const parameter* param){

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

    // If next_edge is still NULL, it means the edge was not found
    if (runner == NULL) {
        next_edge = edgej;
    }
    else{
        next_edge = runner;
    }

    // If prev_edge is still NULL, it means the edge was the first one
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

    const jcv_real Aj = jcv_area(sj);
    const jcv_real Pj = jcv_perimeter(sj);

    jcv_point dEj_dh2 = force_h(Aj, Pj, h7, h2, h3, param);
    jcv_point dEj_dh3 = force_h(Aj, Pj, h2, h3, h8, param);

    jcv_real jacobian_h2[2][2] = {{0, 0}, {0, 0}};
    jcv_real jacobian_h3[2][2] = {{0, 0}, {0, 0}};

    derivative(ri, rl, rj, jacobian_h2);
    derivative(ri, rj, rk, jacobian_h3);

    jcv_real dEj_drix = (dEj_dh2.x*jacobian_h2[0][0] + dEj_dh2.y*jacobian_h2[0][1]
                        + dEj_dh3.x*jacobian_h3[0][0] + dEj_dh3.y*jacobian_h3[0][1]);
    
    jcv_real dEj_driy = (dEj_dh2.x*jacobian_h2[1][0] + dEj_dh2.y*jacobian_h2[1][1]
                        + dEj_dh3.x*jacobian_h3[1][0] + dEj_dh3.y*jacobian_h3[1][1]);

    return (jcv_point){-dEj_drix, -dEj_driy};

}

jcv_point get_edge_force_ii(const jcv_site* si, const parameter* param){
    
    const jcv_point* ri = &(si->p);
    jcv_point* rj = NULL;
    jcv_point* rk = NULL;

    jcv_graphedge* edge = si->edges;
    jcv_graphedge* edge_after = edge->next;


    jcv_real dE_drx = 0;
    jcv_real dE_dry = 0;
    jcv_real jacobian[2][2] = {{0, 0}, {0, 0}};

    jcv_real area = jcv_area(si);
    jcv_real perimeter = jcv_perimeter(si);
    while (edge_after != NULL){

        rj = &(edge->neighbor->p);
        rk = &(edge_after->neighbor->p);

        derivative(ri, rj, rk, jacobian);
        jcv_point dE_dh = force_h(area, perimeter, &(edge->pos[0]), &(edge->pos[1]), &(edge_after->pos[1]), param);
        dE_drx += dE_dh.x*jacobian[0][0] + dE_dh.y*jacobian[0][1];
        dE_dry += dE_dh.x*jacobian[1][0] + dE_dh.y*jacobian[1][1];

        edge = edge_after;
        edge_after = edge_after->next;
    }

    //Last triangle in Delaunay made of last edge and first edge
    derivative(ri, rk, &(si->edges->neighbor->p), jacobian);
    jcv_point dE_dh = force_h(area, perimeter, &(edge->pos[0]), &(edge->pos[1]), &(si->edges->pos[1]), param);
    dE_drx += dE_dh.x*jacobian[0][0] + dE_dh.y*jacobian[0][1];
    dE_dry += dE_dh.x*jacobian[1][0] + dE_dh.y*jacobian[1][1];

    //printf("%f %f \n", si->edges->pos[0].x, edge->pos[1].x);
    return (jcv_point){-dE_drx, -dE_dry};

}
void derivative(const jcv_point* ri, const jcv_point* rj, const jcv_point* rk, jcv_real jacobian[2][2]){
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

    if (D == 0)
       printf("\nri = (%f, %f), rj =  (%f, %f), rk = (%f, %f)\n", ri->x, ri->y, rj->x, rj->y, rk->x, rk->y);

    jcv_real alpha = rjk_sq*jcv_dot(rij, rik)/D;
    jcv_real beta = rik_sq*jcv_dot(rji, rjk)/D;
    jcv_real gamma = rij_sq*jcv_dot(rki, rkj)/D;

    jcv_point dDdri = {4*rjk.y*rij_cross_rjk, -4*rjk.x*rij_cross_rjk}; 

    jcv_point dalphadri = jcv_sub(jcv_mul(rjk_sq/D, jcv_add(rij, rik)), jcv_mul(alpha/D, dDdri));
    jcv_point dbetadri = jcv_sub(jcv_add(jcv_mul(2*jcv_dot(rji, rjk)/D, rik), jcv_mul(-rik_sq/D, rjk)), jcv_mul(beta/D, dDdri));
    jcv_point dgammadri = jcv_sub(jcv_add(jcv_mul(2*jcv_dot(rki, rkj)/D, rij), jcv_mul(-rij_sq/D, rkj)), jcv_mul(gamma/D, dDdri));

    

    // h = alpha*ri + beta*rj + gamma*rk
    // first index = derivative, second index = h
    // ∂h.x/∂ri.x
    jacobian[0][0] = alpha + dalphadri.x*ri->x + dbetadri.x*rj->x + dgammadri.x*rk->x;
    // ∂h.y/∂ri.x  
    jacobian[0][1] = dalphadri.x*ri->y + dbetadri.x*rj->y + dgammadri.x*rk->y;
    // ∂h.x/∂ri.y 
    jacobian[1][0] = dalphadri.y*ri->x + dbetadri.y*rj->x + dgammadri.y*rk->x;
    // ∂h.y/∂ri.y 
    jacobian[1][1] = alpha + dalphadri.y*ri->y + dbetadri.y*rj->y + dgammadri.y*rk->y;
}

void compute_force(data* sys){

    for (int i = 0; i < sys->N_pbc; i++){
        if (sys->sites[i].index >= sys->N) continue;
        sys->forces[sys->sites[i].index] = get_edge_force_ii(sys->sites + i, &sys->parameter);
        jcv_graphedge* graph_edge = sys->sites[i].edges;
        while (graph_edge){
            jcv_point force_ji = get_edge_force_ji(sys->sites + i, graph_edge, &sys->parameter);
            sys->forces[sys->sites[i].index].x += force_ji.x;
            sys->forces[sys->sites[i].index].y += force_ji.y;
            graph_edge = graph_edge->next;
        }
    }
}

jcv_real energy(const jcv_site* sites, const int N, const int N_pbc, const parameter* param){
    jcv_real E = 0;
    for (int i = 0; i < N_pbc; i++){
        if (sites[i].index >= N) continue;
        E += param->Ka*(jcv_area(sites + i) - param->Ao)*(jcv_area(sites + i) - param->Ao);
        E += param->Kp*(jcv_perimeter(sites + i) - param->Po)*(jcv_perimeter(sites + i) - param->Po);
    }
    return E;
}

