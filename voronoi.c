
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define JC_VORONOI_IMPLEMENTATION
// If you wish to use doubles
//#define JCV_REAL_TYPE double
//#define JCV_FABS fabs
//#define JCV_ATAN2 atan2
//#define JCV_CEIL ceilP
//#define JCV_FLOOR floor
//#define JCV_FLT_MAX 1.7976931348623157E+308
#include "jc_voronoi.h"
#include "mersenne.c"

#include "helper.h"
#include "force.h"

int N = 500;
jcv_real L = 1;

int main() {
	 
	//jcv_rect bounding_box = {{-L, L}, {2*L, 2*L}};
	jcv_diagram diagram;
	jcv_point points[9*N];
	const jcv_site* sites;
	jcv_graphedge* graph_edge;

	memset(&diagram, 0, sizeof(jcv_diagram));

	init_genrand(0);
	

	for (int i = 0; i < N; i++) {
		points[i].x = drand(0, L);
		points[i].y = drand(0, L);
	}

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
	
	clock_t start_time, end_time;
    double cpu_time_used;

    start_time = clock();
	jcv_diagram_generate(9*N, points, NULL, 0, &diagram);
	sites = jcv_diagram_get_sites(&diagram);
	end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Took: %f seconds to generate the diagram\n", cpu_time_used);
	
	FILE *file = fopen("edges.txt", "w");
	if (file == NULL) {
		perror("Failed to open file");
		return 1;
	}

	fprintf(file, "%d\n", N);

	for (int i = 0; i < N; i++) {
		fprintf(file, "%f %f\n", (float)points[i].x, (float)points[i].y);
	}

	start_time = clock();
	for (int i = 0; i < 9*N; i++) {
		if (sites[i].index < N){
			graph_edge = sites[i].edges;
			get_edge_force_ii(sites + i);
			while (graph_edge) {
				fprintf(file, "%f %f\n", (float)graph_edge->pos[0].x, (float)graph_edge->pos[0].y);
				fprintf(file, "%f %f\n", (float)graph_edge->pos[1].x, (float)graph_edge->pos[1].y);
				get_edge_force_ji(sites + i, graph_edge);
				graph_edge = graph_edge->next;
			}
		}
	}
	end_time = clock();
	cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	printf("Took: %f seconds to calculate the forces\n", cpu_time_used);

	fclose(file);

	jcv_diagram_free(&diagram);
}
