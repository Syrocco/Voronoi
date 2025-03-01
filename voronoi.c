
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

int N = 30;
int M = 1000;

jcv_real L;
jcv_real dt = 0.01;
FILE *file;

int main(){
	 L = JCV_SQRT((JCV_REAL_TYPE)N);
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

	
	populate_points(points, N, L);
	

	jcv_diagram_generate(9*N, points, NULL, 0, &diagram);
	sites = jcv_diagram_get_sites(&diagram);



	char filename[256];
	jcv_point force[N];

	for (int m = 0; m < M; m++){
		if (m%30 == 0){
			sprintf(filename, "dump/%d.txt", m);
			write(file, filename, points, sites, N);
		}
		for (int i = 0; i < 9*N; i++){
			
			if (sites[i].index >= N) continue;
			force[sites[i].index] = get_edge_force_ii(sites + i);
			graph_edge = sites[i].edges;
			while (graph_edge){
				jcv_point force_ji = get_edge_force_ji(sites + i, graph_edge);
				force[sites[i].index].x += force_ji.x;
				force[sites[i].index].y += force_ji.y;
				graph_edge = graph_edge->next;
			}
		}


		jcv_point f = {0, 0};
		for (int i = 0; i < N; i++){
			points[i].x = pbc(points[i].x + dt*force[i].x, L);
			points[i].y = pbc(points[i].y + dt*force[i].y, L);
			f.x += force[i].x;
			f.y += force[i].y;
		}
		printf("m = %d, f = (%f, %f), E = %f \n", m, f.x, f.y, energy(sites, N));
		populate_points(points, N, L);

		jcv_diagram_generate(9*N, points, NULL, 0, &diagram);
		sites = jcv_diagram_get_sites(&diagram);
	}
}
