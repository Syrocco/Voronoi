#include <stdio.h>
#include "jconfig.h"
#include "helper.h"

/*
Helper functions 
*/

inline jcv_point jcv_add(const jcv_point a, const jcv_point b){
    jcv_point c = {a.x + b.x, a.y + b.y};
    return c;
}

inline jcv_point jcv_sub(const jcv_point a, const jcv_point b){
    jcv_point c = {a.x - b.x, a.y - b.y};
    return c;
}

inline jcv_real jcv_cross(const jcv_point a, const jcv_point b){
    return a.x*b.y - a.y*b.x;
}

inline jcv_real jcv_dot(const jcv_point a, const jcv_point b){
    return a.x*b.x + a.y*b.y;
}

inline jcv_point jcv_mul(const jcv_real a, const jcv_point b){
    jcv_point c = {a*b.x, a*b.y};
    return c;
}

inline jcv_point jcv_vec_mul(const jcv_point a, const jcv_point b){
    jcv_point c = {a.x*b.x, a.y*b.y};
    return c;
}

inline jcv_real jcv_perimeter(const jcv_site* a){
    jcv_real perimeter = 0;
    const jcv_graphedge* e = a->edges;
    while (e){
        perimeter += jcv_point_dist(&e->pos[0], &e->pos[1]);
        e = e->next;
    }
    return perimeter;
}

inline jcv_real jcv_area(const jcv_site* a){
    jcv_real area = 0;
    const jcv_graphedge* e = a->edges;
    while (e){
        area += jcv_cross(e->pos[0], e->pos[1]);
        e = e->next;
    }
    return 0.5*jcv_abs(area);
}

inline jcv_real jcv_lenght(const jcv_point* a){
    return JCV_SQRT(a->x*a->x + a->y*a->y);
}

inline jcv_real jcv_lenght_sq(const jcv_point* a){
    return a->x*a->x + a->y*a->y;
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
    jcv_real rij_sq = jcv_lenght_sq(&rij);
    jcv_real rik_sq = jcv_lenght_sq(&rik);
    jcv_real D = 2*rij_cross_rjk*rij_cross_rjk;


    jcv_real alpha = rjk_sq*jcv_dot(rij, rik)/D;
    jcv_real beta = rik_sq*jcv_dot(rji, rjk)/D;
    jcv_real gamma = rij_sq*jcv_dot(rki, rkj)/D;

    //Derivative of D
    jcv_real point = jcv_cross(rij, rjk);
    jcv_point dDdri = {4*rjk.y*point, -4*rjk.x*point};

    jcv_point dalphadri = jcv_sub(jcv_mul(rjk_sq/D, jcv_add(rij, rik)), jcv_mul(alpha/D, dDdri));
    jcv_point dbetadri = jcv_sub(jcv_add(jcv_mul(2*jcv_dot(rji, rjk)/D, rik), jcv_mul(-rik_sq/D, rjk)), jcv_mul(beta/D, dDdri));
    jcv_point dgammadri = jcv_sub(jcv_add(jcv_mul(2*jcv_dot(rki, rkj)/D, rij), jcv_mul(-rij_sq/D, rkj)), jcv_mul(gamma/D, dDdri));


    //first index = derivative, second index = h
    jacobian[0][0] = alpha + dalphadri.x*ri->x + dbetadri.x*rj->x + dgammadri.x*rk->x;
    jacobian[0][1] = dalphadri.x*ri->y + dbetadri.x*rj->y + dgammadri.x*rk->y;
    jacobian[1][0] = dalphadri.y*ri->x + dbetadri.y*rj->x + dgammadri.y*rk->x;
    jacobian[1][1] = alpha + dalphadri.y*ri->y + dbetadri.y*rj->y + dgammadri.y*rk->y;
}