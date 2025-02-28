#include "force.h"
#include "helper.h"
#include "jc_voronoi.h"
#include <stdio.h>

jcv_real Ao = 1.0;
jcv_real Po = 1.0;
jcv_real Ka = 1.0;
jcv_real Kp = 1.0;


jcv_point force_h(const jcv_real A, const jcv_real P, const jcv_point* h7, const jcv_point* h2, const jcv_point* h3){
    jcv_real h72 = jcv_point_dist(h2, h7);
    jcv_real h23 = jcv_point_dist(h2, h3);
    jcv_point opp = {h3->y - h7->y, h3->x - h7->x};
    jcv_point area = jcv_mul(Ka*(A - Ao), opp);

    jcv_point perimeter_last = jcv_add(jcv_mul(1/h72, jcv_sub(*h2, *h7)), jcv_mul(1/h23, jcv_sub(*h2, *h3)));

    jcv_point perimeter = jcv_mul(2*Kp*(P - Po), perimeter_last);
    return jcv_add(area, perimeter);
}

jcv_point get_edge_force_ji(const jcv_site* si, const jcv_graphedge* edgei){

    //edgei is owned by site i: si

    jcv_site* sj = edgei->neighbor;
    jcv_graphedge* edgej = sj->edges;

    

    jcv_graphedge* runner = edgej;

    jcv_graphedge* prev_edge = NULL;
    jcv_graphedge* next_edge = NULL;

    while (runner->next != NULL){
        if (runner->next->neighbor == si) {
            runner = runner->next;
            break;
        }
        prev_edge = runner->next;
        runner = runner->next;
    }

    //If runner not at the end, h3-h8 is runner->next
    if (runner->next != NULL){
        next_edge = runner->next;
    }
    // If at the end, h3-h8 was the initial starting point
    else{
        next_edge = edgej;
    }

    //If prev_edge is NULL, h2-h7 is the last edge
    if (prev_edge == NULL){
        while (runner->next != NULL){
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
    const jcv_point* rk = &(next_edge->neighbor->p);
    const jcv_point* rl = &(prev_edge->neighbor->p);

    const jcv_real Aj = jcv_area(sj);
    const jcv_real Pj = jcv_perimeter(sj);

    
    
    //printf("Sites: ri: %f %f\n rj: %f %f\n rk: %f %f\n rl: %f %f\n\n", ri->x, ri->y, rj->x, rj->y, rk->x, rk->y, rl->x, rl->y);
    //printf("Edges: h2: %f %f\n h3: %f %f\n h7: %f %f\n h8: %f %f\n\n", h2->x, h2->y, h3->x, h3->y, h7->x, h7->y, h8->x, h8->y);




    jcv_real Fj = 0;
    jcv_real h72 = jcv_point_dist(h2, h7);
    jcv_real h23 = jcv_point_dist(h2, h3);
    jcv_real h83 = jcv_point_dist(h3, h8);
    jcv_point dEj_dh2 = force_h(Aj, Pj, h7, h2, h3);
    jcv_point dEj_dh3 = force_h(Aj, Pj, h2, h3, h8);

    jcv_real jacobian_h2[2][2] = {{0, 0}, {0, 0}};
    jcv_real jacobian_h3[2][2] = {{0, 0}, {0, 0}};
    derivative(ri, rl, rj, jacobian_h2);
    derivative(ri, rj, rk, jacobian_h3);

    jcv_real dEj_drix = (dEj_dh2.x*jacobian_h2[0][0] + dEj_dh2.y*jacobian_h2[0][1]
                        + dEj_dh3.x*jacobian_h3[0][0] + dEj_dh3.y*jacobian_h3[0][1]);
    
    jcv_real dEj_driy = (dEj_dh2.x*jacobian_h2[1][0] + dEj_dh2.y*jacobian_h2[1][1]
                        + dEj_dh3.x*jacobian_h3[1][0] + dEj_dh3.y*jacobian_h3[1][1]);

    return (jcv_point){dEj_drix, dEj_driy};

}

jcv_point get_edge_force_ii(const jcv_site* si){
    
    const jcv_point* ri = &(si->p);
    jcv_point* rj,* rk;

    jcv_graphedge* edge = si->edges;
    jcv_graphedge* edge_after = edge->next;


    jcv_real dE_drx = 0;
    jcv_real dE_dry = 0;
    jcv_real jacobian[2][2] = {{0, 0}, {0, 0}};

    while (edge_after != NULL){

        rj = &(edge->neighbor->p);
        rk = &(edge_after->neighbor->p);

        derivative(ri, rj, rk, jacobian);
        jcv_point dE_dh = force_h(jcv_area(si), jcv_perimeter(si), &(edge->pos[0]), &(edge->pos[1]), &(edge_after->pos[1]));
        dE_drx += dE_dh.x*jacobian[0][0] + dE_dh.y*jacobian[0][1];
        dE_drx += dE_dh.x*jacobian[1][0] + dE_dh.y*jacobian[1][1];

        edge = edge_after;
        edge_after = edge_after->next;
    }

    //Last triangle in Delaunay made of last edge and first edge
    derivative(ri, rk, &(si->edges->neighbor->p), jacobian);
    jcv_point dE_dh = force_h(jcv_area(si), jcv_perimeter(si), &(edge->pos[0]), &(edge->pos[1]), &(si->edges->pos[1]));
    dE_drx += dE_dh.x*jacobian[0][0] + dE_dh.y*jacobian[0][1];
    dE_drx += dE_dh.x*jacobian[1][0] + dE_dh.y*jacobian[1][1];

    //printf("%f %f \n", si->edges->pos[0].x, edge->pos[1].x);
    return (jcv_point){dE_drx, dE_dry};

}