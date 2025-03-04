#include "voronoi.h"
#include "jc_voronoi.h"
#include "helper.h"
#include <math.h>

void randomInitial(data* sys){
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++){
        sys->positions[i].x = drand(0, sys->L);
        sys->positions[i].y = drand(0, sys->L);
        addBoundary(sys, i);
    }
}

void rsaInitial(data* sys, jcv_real min_distance) {
    int max_attempts = 1000;
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        int attempts = 0;
        int placed = 0;
        while (attempts < max_attempts && !placed) {
            double x = drand(0, sys->L);
            double y = drand(0, sys->L);
            int overlap = 0;
            for (int j = 0; j < sys->N; j++) {
                double dx = x - sys->positions[j].x;
                double dy = y - sys->positions[j].y;
                if (sqrt(dx * dx + dy * dy) < min_distance) {
                    overlap = 1;
                    break;
                }
            }
            if (!overlap) {
                sys->positions[i].x = x;
                sys->positions[i].y = y;
                addBoundary(sys, i);
                placed = 1;
            }
            attempts++;
        }
        if (!placed) {
            printf("Failed to place particle %d after %d attempts\n", i, max_attempts);
        }
    }
}