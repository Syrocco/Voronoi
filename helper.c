#include <stdio.h>
#include "jconfig.h"
#include "helper.h"
#include <stdio.h>
#include "voronoi.h"
/*
Helper functions 
*/
#define eps 0.000001

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


void populate_points(jcv_point* points, int N, jcv_real L){
    int count = 0;
	for (int i = 0; i < N; i++) {
		for (int j = -1; j < 2; j++){
			for (int k = -1; k < 2; k++){
				if (j == 0 && k == 0) continue;
				points[N + count].x = points[i].x + j*L;
				points[N + count].y = points[i].y + k*L;
				count++;
			}
		}
	}
}

void pbc(jcv_point* p, jcv_real L, jcv_real gamma){
    if (p->y < 0){
        p->y += L;
        p->x += gamma*L;
    }
    else if (p->y > L){
        p->y -= L;
        p->x -= gamma*L;
    }
    
    if (p->x - p->y*gamma < 0){
        p->x += L;
    }
    else if (p->x - p->y*gamma > L){
        p->x -= L;
    }
}

//Will not work with dL > 0.5
void addBoundary(data* sys, int i){
    /* int val = 1;
    for (int l = -val; l <= val; l++){
        for (int j = -val; j <= val; j++){
            if (l == 0 && j == 0) continue;
            sys->positions[sys->N_pbc].x = sys->positions[i].x + l*sys->L;
            sys->positions[sys->N_pbc].y = sys->positions[i].y + j*sys->L;
            sys->N_pbc++;
        }
    } */

    jcv_real gamma = sys->gamma;
    jcv_real L = sys->L;
    jcv_real dL = sys->dL;

    // Adjust the conditions for left and right boundaries
    int left = (sys->positions[i].x - sys->positions[i].y*gamma)/L < dL;
    int right = (sys->positions[i].x - sys->positions[i].y*gamma)/L > 1 - dL;
    int bottom = sys->positions[i].y/L < dL;
    int top = sys->positions[i].y/L >  1 - dL;
    if (left || right) {
        sys->positions[sys->N_pbc].x = sys->positions[i].x + (left ? L : -L);
        sys->positions[sys->N_pbc].y = sys->positions[i].y;
        sys->N_pbc++;
    } 
    if (bottom || top) {
        sys->positions[sys->N_pbc].x = sys->positions[i].x + gamma*(bottom ? L : -L);
        sys->positions[sys->N_pbc].y = sys->positions[i].y + (bottom ? L : -L);
        sys->N_pbc++;
    } 
    if ((left && bottom) || (right && bottom) || (left && top) || (right && top)) {
        sys->positions[sys->N_pbc].x = sys->positions[i].x + gamma*(bottom ? L : -L) + (left ? L : -L);
        sys->positions[sys->N_pbc].y = sys->positions[i].y + (bottom ? L : -L);
        sys->N_pbc++;
    } 
}

//Box Muller gaussian number generator
void gaussian(jcv_real* n1, jcv_real* n2){
    jcv_real u1 = drand(0, 1);
    jcv_real u2 = drand(0, 1);
    *n1 = JCV_SQRT(-2*JCV_LOG(u1))*JCV_COS(2*JCV_PI*u2);
    *n2 = JCV_SQRT(-2*JCV_LOG(u1))*JCV_SIN(2*JCV_PI*u2);
}