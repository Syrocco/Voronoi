#include "voronoi.h"
#include "jc_voronoi.h"
#include "helper.h"
#include <math.h>
#include "parser.h"
#include "logger.h"

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

void distribute_area(data* sys){
    for (int i = 0; i < sys->N; i++){
        jcv_real fac = (((jcv_real) i)/sys->N < sys->n_frac_small ? 1 : sys->size_large_over_small);
        sys->prefered_area[i] = fac/(sys->n_frac_small + sys->size_large_over_small - sys->n_frac_small*sys->size_large_over_small);
    }
}

void read_from_dump_initial(data* sys, char* filename, int frame){
    
    Dump* dump = dump_open(filename,'r');		
    sys->N = get_natoms(dump);

    if (frame < 0){
        frame = dump->nframes + frame;
    }
    jump_to_frame(frame, dump);

    char hboxy, hboxx;
    sys->L = get_boxy(&hboxy, 1, dump);
    sys->gamma = (get_boxx(&hboxx, 0, dump) - sys->L)/sys->L;
    sys->parameter.gamma_rate = 0;



    double x[sys->N], y[sys->N], prefered_area[sys->N];
    
    get_doubleatomprop("x", x, sys->N, dump);
    get_doubleatomprop("y", y, sys->N, dump);
    char header[LINESIZE];
    get_header("ATOMS", header, LINESIZE, dump);
    //check if header contains prefered_area
    if (strstr(header, "prefered_area") != NULL){
        get_doubleatomprop("prefered_area", prefered_area, sys->N, dump);
    }
    else{
        distribute_area(sys);
    }

    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++){
        sys->positions[i].x = x[i];
        sys->positions[i].y = y[i];
        sys->prefered_area[i] = prefered_area[i];
        addBoundary(sys, i);
    }

    
    dump_close(dump);

}