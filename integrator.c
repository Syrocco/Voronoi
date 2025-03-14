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


jcv_real line_search(data* sys, jcv_point* gradient, jcv_point *direction, double alpha, double rho, double c){
    jcv_point x_old[sys->N];

    double dot_product = 0.0;

    jcv_real old_energy = energy_total(sys);

    // Compute dot product of gradient and search direction
    for (int i = 0; i < sys->N; i++) {
        x_old[i] = sys->positions[i];
        dot_product += jcv_dot(gradient[i], direction[i]);
    }
    //dot_product /= sys->N;
    if (dot_product > 0){
        dot_product = 0;
        for (int i = 0; i < sys->N; i++) {
            direction[i] = sys->forces[i];
            dot_product += jcv_dot(gradient[i], direction[i]);
        }
    }

    while (1) {
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = x_old[i].x + alpha*direction[i].x;
            sys->positions[i].y = x_old[i].y + alpha*direction[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }

        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        jcv_real new_energy = energy_total(sys);
        //printf("new_energy = %.9lf, old_energy = %.9lf, alpha = %.9lf\n", new_energy, old_energy, alpha);

        if (new_energy <= old_energy + c * alpha * dot_product) {
            break;  // Stop if condition is satisfied
        }
        alpha *= rho;  // Reduce step size
    }
    return alpha;
}

void conjugateGradientStep(data* sys) {
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);

    compute_force(sys);

    loggers(sys);

    jcv_real chi = 0.0;
    jcv_point gradient[sys->N];
    jcv_point delta[sys->N];

    for (int i = 0; i < sys->N; i++) {
        gradient[i].x = -sys->forces[i].x; // grad(U)
        gradient[i].y = -sys->forces[i].y;
        delta[i] = sys->forces[i]; // -grad(U)
    }

    

    int count = 0;
    while (1) {
        count++;
        line_search(sys, gradient, delta, 1, 0.5, 0.0001);

        compute_force(sys);
        jcv_real gnorm = 0.0;
        for (int i = 0; i < 2 * sys->N; i++) {
            gnorm += jcv_lenght_sq(&sys->forces[i]);
        }
        gnorm = JCV_SQRT(gnorm/sys->N);
        //printf("gnorm = %.14lf\n", gnorm);
        if (gnorm < 1e-6) {
            break;
        }
        else{
            jcv_real numer = 0;
            jcv_real denom = 0;
            for (int i = 0; i < sys->N; i++) {
                numer += jcv_dot(jcv_add(sys->forces[i], gradient[i]), sys->forces[i]);
                denom += jcv_lenght_sq(&gradient[i]);
            }
            chi = jcv_max(0., numer/denom);
            for (int i = 0; i < sys->N; i++) {
                gradient[i].x = -sys->forces[i].x;
                gradient[i].y = -sys->forces[i].y;
                delta[i].x = sys->forces[i].x + chi*delta[i].x;
                delta[i].y = sys->forces[i].y + chi*delta[i].y;
            }
        }
        //saveTXT(sys);
    }
    shear(sys);
    if (sys->i - sys->shear_start >= 0) {
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = sys->positions[i].x + sys->parameter.gamma_rate * sys->dt * sys->positions[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
    }
} 

void fireStep(data* sys){
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);

    loggers(sys);



    jcv_real P = 0.0;
    jcv_real fnorm = 0.0;
    jcv_real vnorm = 0.0;
    
    
    const jcv_real alpha_start = 0.3;
   
    const jcv_real f_inc = 1.1;
    const jcv_real f_dec = 0.5;
    const jcv_real f_alpha = 0.99;
    const int n_delay = 20;
    
    const jcv_real dt_max = sys->dt_fire;
    const jcv_real dt_min = sys->dt_fire/100;

    int n_pos = 0;
    jcv_real dt_now = sys->dt_fire;
    jcv_real alpha_now = alpha_start;

    for (int i = 0; i < sys->N; i++){
        sys->velocities[i].x = 0;
        sys->velocities[i].y = 0;
    }
    int count = 0;
    do{
        

        count++;
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
        if (count > 2500){
            if (fnorm > 0.01)
                dt_now = dt_max*1.5;
            printf("fnorm = %.14lf, vnorm = %lf, P = %lf, dt = %lf, E = %lf\n", fnorm, vnorm, P, dt_now, energy_total(sys));
            //saveTXT(sys);
            if (count > 20000){
                printf("Failed to converge\n");
                unexpectedClosure(sys);
                exit(3);
            } 
        }


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
  

        
    } while (fnorm/sys->L > 1e-8); // L = sqrt((float)N)
    shear(sys); 
    if (sys->i - sys->shear_start >= 0){
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++){
            sys->positions[i].x = sys->positions[i].x + sys->parameter.gamma_rate*sys->dt*sys->positions[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
    }

}

