#include "integrator.h"

#include "helper.h"
#include "force.h"
#include "voronoi.h"
#include "logger.h"
#include "thermo.h"
#include <time.h>
#include <omp.h>

#define TIME_FUNCTION(func, ...) \
    do { \
        clock_t start_time = clock(); \
        func(__VA_ARGS__); \
        clock_t end_time = clock(); \
        double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC; \
        printf("Time for %s: %f seconds\n", #func, elapsed_time); \
    } while (0)




void eulerStep(data* sys){


    //TIME_FUNCTION(jcv_diagram_generate,sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    


    loggers(sys);

        
    shear(sys);

    //TIME_FUNCTION(compute_force, sys);
    compute_force(sys);

    sys->N_pbc = sys->N;
    jcv_real forcex = 0;
    jcv_real forcey = 0;
    for (int i = 0; i < sys->N; i++){
        #if NOISE 
        jcv_real r1, r2;
        gaussian(&r1, &r2);
        sys->positions[i].x = sys->positions[i].x + sys->dt*sys->forces[i].x + JCV_SQRT(2*sys->parameter.T*sys->dt)*r1; 
        sys->positions[i].y = sys->positions[i].y + sys->dt*sys->forces[i].y + JCV_SQRT(2*sys->parameter.T*sys->dt)*r2;
        #else
        sys->positions[i].x = sys->positions[i].x + sys->dt*sys->forces[i].x; 
        sys->positions[i].y = sys->positions[i].y + sys->dt*sys->forces[i].y;
        #endif
        if (sys->i - sys->shear_start >= 0){
            sys->positions[i].x = sys->positions[i].x + sys->parameter.gamma_rate*sys->dt*sys->positions[i].y;
        }
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
        forcex += sys->forces[i].x;
        forcey += sys->forces[i].y;
    }
    
    if ((forcex*forcex + forcey*forcey)/sys->N > 0.00001){
        //printf("forcex = %.9lf, forcey = %.9lf \n", forcex, forcey);
        //exit(3);
    }

}

void fireStep(data* sys){
    if ((sys->i + 1)%100 == 0){
        sys->parameter.gamma_rate *= -1;
    }
    sys->gamma += sys->dt*sys->parameter.gamma_rate;
    jcv_real P = 0.0;
    jcv_real fnorm = 0.0;
    jcv_real vnorm = 0.0;
    
    
    const jcv_real alpha_start = 0.3;
   
    const jcv_real f_inc = 1.1;
    const jcv_real f_dec = 0.5;
    const jcv_real f_alpha = 0.99;
    const int n_delay = 5;
    
    const jcv_real dt_max = sys->dt;
    const jcv_real dt_min = sys->dt/100;

    int n_pos = 0;
    jcv_real dt_now = sys->dt;
    jcv_real alpha_now = alpha_start;

    for (int i = 0; i < sys->N; i++){
        sys->velocities[i].x = 0;
        sys->velocities[i].y = 0;
    }
    do{

        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, 0, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);

        compute_force(sys);

        P = 0;
        fnorm = 0;
        vnorm = 0;

        for (int i = 0; i < sys->N; i++) {
            sys->velocities[i].x += dt_now*sys->forces[i].x;
            sys->velocities[i].y += dt_now*sys->forces[i].y;
            P += sys->forces[i].x*sys->velocities[i].x + sys->forces[i].y*sys->velocities[i].y;
            fnorm += sys->forces[i].x*sys->forces[i].x + sys->forces[i].y*sys->forces[i].y;
            vnorm += sys->velocities[i].x*sys->velocities[i].x + sys->velocities[i].y*sys->velocities[i].y;
        }

        fnorm = JCV_SQRT(fnorm);
        vnorm = JCV_SQRT(vnorm);

        printf("fnorm = %.10lf, dt_now = %f, P = %.10lf, alpha_now = %f \n", fnorm/sys->L, dt_now, P, alpha_now);


        if (P > 0){
            n_pos++;
            for (int i = 0; i < sys->N; i++) {
                sys->velocities[i].x = (1.0 - alpha_now)*sys->velocities[i].x + alpha_now*sys->forces[i].x*(vnorm/fnorm);
                sys->velocities[i].y = (1.0 - alpha_now)*sys->velocities[i].y + alpha_now*sys->forces[i].y*(vnorm/fnorm);
            }
            if (n_pos > n_delay) {
                dt_now = jcv_min(dt_now*f_inc, dt_max);
                alpha_now *= f_alpha;
            }
        } 
        else{
            n_pos = 0;
            dt_now = jcv_max(dt_now*f_dec, dt_min);
            alpha_now = alpha_start;
            
            // Zero the velocities
            for (int i = 0; i < sys->N; i++){
                sys->velocities[i].x = 0;
                sys->velocities[i].y = 0;
            }
        }
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x += dt_now*sys->velocities[i].x;
            sys->positions[i].y += dt_now*sys->velocities[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        
        loggers(sys);
        if (fnorm/sys->L < 1e-9){
            printf("%lf \n", (fnorm/sys->L)/1e-10);
        }
    } while (fnorm/sys->L > 1e-10); // L = sqrt((float)N)
   
    exit(3);
}

