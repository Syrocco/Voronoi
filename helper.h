#ifndef HELPER_H
#define HELPER_H

#include "jc_voronoi.h"

jcv_real jcv_perimeter(const jcv_site* a);
jcv_real jcv_area(const jcv_site* a);

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
        area += e->pos[0].x*e->pos[1].y - e->pos[0].y*e->pos[1].x;
        e = e->next;
    }
    return 0.5*jcv_abs(area);
}

#endif // HELPER_H