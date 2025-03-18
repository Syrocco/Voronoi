#include "integrator.h"

#include "helper.h"
#include "force.h"
#include "voronoi.h"
#include "logger.h"
#include "thermo.h"
#include <time.h>
#include <omp.h>

#define STAY_IN_LOCAL_MINIMA 0
#define FORCE_BASED_MINIMIZATION 1

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

    //TIME_FUNCTION(compute_force, sys);
    compute_force(sys);
    
    //Shear applied after the force.
    shear(sys);
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
    sys->time += sys->dt;
    
    if ((forcex*forcex + forcey*forcey)/sys->N > 0.00001){
        //printf("forcex = %.9lf, forcey = %.9lf \n", forcex, forcey);
        //exit(3);
    }
}

void backwardEulerStep(data* sys) {
    jcv_point initial_positions[sys->N];
    jcv_real tolerance = 1e-6;
    jcv_real error = 1.0;
    int max_iterations = 50;
    int iter = 0;
    
    // Generate initial Voronoi diagram and log data
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    loggers(sys);
    
    // Store initial positions
    for (int i = 0; i < sys->N; i++) {
        initial_positions[i] = sys->positions[i];
    }
    
    // Initial guess using explicit Euler
    compute_force(sys);
    
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        sys->positions[i].x = initial_positions[i].x + sys->dt * sys->forces[i].x;
        sys->positions[i].y = initial_positions[i].y + sys->dt * sys->forces[i].y;
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
    }
    
    // Newton iterations to solve the implicit equation
    while (error > tolerance && iter < max_iterations) {
        // Compute forces at current position estimate
        
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        
        // Compute residual and update positions
        error = 0.0;
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            // Residual: r = x - x_0 - h*f(x)
            jcv_real rx = sys->positions[i].x - initial_positions[i].x - sys->dt * sys->forces[i].x;
            jcv_real ry = sys->positions[i].y - initial_positions[i].y - sys->dt * sys->forces[i].y;
            
            // Update position (x_new = x - r)
            sys->positions[i].x -= rx;
            sys->positions[i].y -= ry;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
            
            // Accumulate error
            error += rx*rx + ry*ry;
        }
        error = sqrt(error / sys->N);
        iter++;
    }
    
    // Apply shear if needed
    shear(sys);
    sys->N_pbc = sys->N;
    if (sys->i - sys->shear_start >= 0) {
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = sys->positions[i].x + sys->parameter.gamma_rate * sys->dt * sys->positions[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
    }
    
    sys->time += sys->dt;
}

void rk4Step(data* sys) {
    jcv_point k1[sys->N], k2[sys->N], k3[sys->N], k4[sys->N];
    jcv_point initial_positions[sys->N];

    for (int i = 0; i < sys->N; i++) {
        initial_positions[i] = sys->positions[i];
    }

    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    loggers(sys);


    // Compute k1 ------------------------

    compute_force(sys);
    for (int i = 0; i < sys->N; i++) {
        k1[i].x = sys->dt * sys->forces[i].x;
        k1[i].y = sys->dt * sys->forces[i].y;
    }
    // ###################################

    // Compute k2 ------------------------
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        sys->positions[i].x = initial_positions[i].x + 0.5 * k1[i].x;
        sys->positions[i].y = initial_positions[i].y + 0.5 * k1[i].y;
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
    }
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    compute_force(sys);
    for (int i = 0; i < sys->N; i++) {
        k2[i].x = sys->dt * sys->forces[i].x;
        k2[i].y = sys->dt * sys->forces[i].y;
    }
    // ###################################

    // Compute k3 ------------------------
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        sys->positions[i].x = initial_positions[i].x + 0.5 * k2[i].x;
        sys->positions[i].y = initial_positions[i].y + 0.5 * k2[i].y;
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
    }
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    compute_force(sys);
    for (int i = 0; i < sys->N; i++) {
        k3[i].x = sys->dt * sys->forces[i].x;
        k3[i].y = sys->dt * sys->forces[i].y;
    }
    // ###################################

    // Compute k4 ------------------------
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        sys->positions[i].x = initial_positions[i].x + k3[i].x;
        sys->positions[i].y = initial_positions[i].y + k3[i].y;
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
    }
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    compute_force(sys);
    for (int i = 0; i < sys->N; i++) {
        k4[i].x = sys->dt * sys->forces[i].x;
        k4[i].y = sys->dt * sys->forces[i].y;
    }
    // ###################################

    // Update positions
    shear(sys);
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        sys->positions[i].x = initial_positions[i].x + (k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x) / 6.0;
        sys->positions[i].y = initial_positions[i].y + (k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y) / 6.0;
        if (sys->i - sys->shear_start >= 0) {
            sys->positions[i].x = sys->positions[i].x + sys->parameter.gamma_rate * sys->dt * sys->positions[i].y;
        }
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
    } 
    sys->time += sys->dt;
}

void rkf45Step(data* sys) {
    jcv_point k1[sys->N], k2[sys->N], k3[sys->N], k4[sys->N], k5[sys->N], k6[sys->N];
    jcv_point initial_positions[sys->N];
    jcv_real error, dt_new;
    const jcv_real tolerance = 1e-4;
    const jcv_real safety = 0.9;
    const jcv_real min_dt = 1e-10;
    const jcv_real max_dt = 1.0;

    for (int i = 0; i < sys->N; i++) {
        initial_positions[i] = sys->positions[i];
    }

    
    jcv_real old_dt = sys->dt;
    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);
    loggers(sys);

    // Compute k1
    compute_force(sys);
    for (int i = 0; i < sys->N; i++) {
        k1[i].x = sys->dt * sys->forces[i].x;
        k1[i].y = sys->dt * sys->forces[i].y;
    }

    int count = 0;
    do {
        
        if (count > 0){
            for (int i = 0; i < sys->N; i++){
                k1[i].x = k1[i].x/old_dt*sys->dt;
                k1[i].y = k1[i].y/old_dt*sys->dt;
            }
        }

        count++;
        // Compute k2
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = initial_positions[i].x + 0.25 * k1[i].x;
            sys->positions[i].y = initial_positions[i].y + 0.25 * k1[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        for (int i = 0; i < sys->N; i++) {
            k2[i].x = sys->dt * sys->forces[i].x;
            k2[i].y = sys->dt * sys->forces[i].y;
        }

        // Compute k3
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = initial_positions[i].x + (3.0 / 32.0) * k1[i].x + (9.0 / 32.0) * k2[i].x;
            sys->positions[i].y = initial_positions[i].y + (3.0 / 32.0) * k1[i].y + (9.0 / 32.0) * k2[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        for (int i = 0; i < sys->N; i++) {
            k3[i].x = sys->dt * sys->forces[i].x;
            k3[i].y = sys->dt * sys->forces[i].y;
        }

        // Compute k4
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = initial_positions[i].x + (1932.0 / 2197.0) * k1[i].x - (7200.0 / 2197.0) * k2[i].x + (7296.0 / 2197.0) * k3[i].x;
            sys->positions[i].y = initial_positions[i].y + (1932.0 / 2197.0) * k1[i].y - (7200.0 / 2197.0) * k2[i].y + (7296.0 / 2197.0) * k3[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        for (int i = 0; i < sys->N; i++) {
            k4[i].x = sys->dt * sys->forces[i].x;
            k4[i].y = sys->dt * sys->forces[i].y;
        }

        // Compute k5
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = initial_positions[i].x + (439.0 / 216.0) * k1[i].x - 8.0 * k2[i].x + (3680.0 / 513.0) * k3[i].x - (845.0 / 4104.0) * k4[i].x;
            sys->positions[i].y = initial_positions[i].y + (439.0 / 216.0) * k1[i].y - 8.0 * k2[i].y + (3680.0 / 513.0) * k3[i].y - (845.0 / 4104.0) * k4[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        for (int i = 0; i < sys->N; i++) {
            k5[i].x = sys->dt * sys->forces[i].x;
            k5[i].y = sys->dt * sys->forces[i].y;
        }

        // Compute k6
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = initial_positions[i].x - (8.0 / 27.0) * k1[i].x + 2.0 * k2[i].x - (3544.0 / 2565.0) * k3[i].x + (1859.0 / 4104.0) * k4[i].x - (11.0 / 40.0) * k5[i].x;
            sys->positions[i].y = initial_positions[i].y - (8.0 / 27.0) * k1[i].y + 2.0 * k2[i].y - (3544.0 / 2565.0) * k3[i].y + (1859.0 / 4104.0) * k4[i].y - (11.0 / 40.0) * k5[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        for (int i = 0; i < sys->N; i++) {
            k6[i].x = sys->dt * sys->forces[i].x;
            k6[i].y = sys->dt * sys->forces[i].y;
        }

        // Estimate error
        error = 0.0;
        for (int i = 0; i < sys->N; i++) {
            jcv_real err_x = (1.0 / 360.0) * k1[i].x - (128.0 / 4275.0) * k3[i].x - (2197.0 / 75240.0) * k4[i].x + (1.0 / 50.0) * k5[i].x + (2.0 / 55.0) * k6[i].x;
            jcv_real err_y = (1.0 / 360.0) * k1[i].y - (128.0 / 4275.0) * k3[i].y - (2197.0 / 75240.0) * k4[i].y + (1.0 / 50.0) * k5[i].y + (2.0 / 55.0) * k6[i].y;
            error += err_x * err_x + err_y * err_y;
        }
        error = sqrt(error / sys->N);

        // Adjust timestep
        dt_new = safety * sys->dt * pow(tolerance / error, 0.2);
        if (dt_new < min_dt) dt_new = min_dt;
        if (dt_new > max_dt) dt_new = max_dt;
        sys->dt = dt_new;

    } while (error > tolerance);


    // Update positions
    shear(sys);
    sys->N_pbc = sys->N;
    for (int i = 0; i < sys->N; i++) {
        sys->positions[i].x = initial_positions[i].x + (16.0 / 135.0) * k1[i].x + (6656.0 / 12825.0) * k3[i].x + (28561.0 / 56430.0) * k4[i].x - (9.0 / 50.0) * k5[i].x + (2.0 / 55.0) * k6[i].x;
        sys->positions[i].y = initial_positions[i].y + (16.0 / 135.0) * k1[i].y + (6656.0 / 12825.0) * k3[i].y + (28561.0 / 56430.0) * k4[i].y - (9.0 / 50.0) * k5[i].y + (2.0 / 55.0) * k6[i].y;
        if (sys->i - sys->shear_start >= 0) {
            sys->positions[i].x = sys->positions[i].x + sys->parameter.gamma_rate * sys->dt * sys->positions[i].y;
        }
        pbc(&sys->positions[i], sys->L, sys->gamma);
        addBoundary(sys, i);
    }

    sys->time += sys->dt;
}

// Line search with backtracking and goldstein-armijo condition
jcv_real line_search(data* sys, jcv_point* gradient, jcv_point *direction, jcv_real alpha, jcv_real rho, jcv_real c){
    jcv_point x_old[sys->N];

    jcv_real dot_product = 0.0;

    jcv_real old_energy = energy_total(sys);

    // Compute dot product of gradient and search direction
    for (int i = 0; i < sys->N; i++) {
        x_old[i] = sys->positions[i];
        dot_product += jcv_dot(gradient[i], direction[i]);
    }

    // If going uphill, change direction to -gradient (pure gradient descent)
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
            break;  // Stop if a Armijo condition is satisfied
        }
        alpha *= rho;  // Reduce step size
    }
    return alpha;
}

// Line search based on force magnitude reduction
jcv_real line_search_by_force(data* sys, jcv_point* gradient, jcv_point* direction, 
    jcv_real alpha, jcv_real rho, jcv_real c __attribute__((unused))) {


    jcv_point x_old[sys->N];

    // Store initial positions
    for (int i = 0; i < sys->N; i++) {
        x_old[i] = sys->positions[i];
    }

    // Calculate initial force magnitude squared
    jcv_real initial_force_mag_sq = 0.0;
    for (int i = 0; i < sys->N; i++) {
        initial_force_mag_sq += jcv_lenght_sq(&sys->forces[i]);
    }

    // Calculate initial dot product to ensure we're going downhill
    jcv_real dot_product = 0.0;
    for (int i = 0; i < sys->N; i++) {
        dot_product += jcv_dot(gradient[i], direction[i]);
    }

    // If going uphill, change direction to -gradient (pure gradient descent)
    if (dot_product > 0) {
        dot_product = 0;
        for (int i = 0; i < sys->N; i++) {
            direction[i] = sys->forces[i];
            dot_product += jcv_dot(gradient[i], direction[i]);
        }
    }

    int max_iterations = 20;  // Prevent infinite loop
    int iter = 0;

    while (iter++ < max_iterations) {
        // Take a step
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = x_old[i].x + alpha * direction[i].x;
            sys->positions[i].y = x_old[i].y + alpha * direction[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }

        // Compute new forces
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);

        // Calculate new force magnitude
        jcv_real new_force_mag_sq = 0.0;
        for (int i = 0; i < sys->N; i++) {
            new_force_mag_sq += jcv_lenght_sq(&sys->forces[i]);
        }

        // Check if force magnitude has decreased
        if (new_force_mag_sq < initial_force_mag_sq) {
            return alpha;  // Accept this step size
        }

        alpha *= rho;  // Reduce step size

    }

    return alpha;  // Return the last tried alpha
}



// Custom mine search with backtracking and goldstein-armijo condition, as well as trying to say in the local minima.
jcv_real line_search_local_minima(data* sys, jcv_point* gradient, jcv_point *direction, jcv_real alpha, jcv_real rho, jcv_real c){
    int can_increase = 1;
    jcv_point x_old[sys->N];

    jcv_real dot_product = 0.0;

    jcv_real old_energy = energy_total(sys);

    // Compute dot product of gradient and search direction
    for (int i = 0; i < sys->N; i++) {
        x_old[i] = sys->positions[i];
        dot_product += jcv_dot(gradient[i], direction[i]);
    }

    // If going uphill, change direction to -gradient (pure gradient descent)
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
        //printf("new_energy = %.16lf, old_energy = %.16lf, alpha = %.9lf\n", new_energy, old_energy, alpha);
        int validity_test = new_energy <= old_energy + c*alpha*dot_product;
        if (validity_test && can_increase == 0) {
            break;  // Stop if a Armijo condition is satisfied
        }
        else if (validity_test && can_increase == 1){
            alpha /= rho; //Start with small alpha, increases it while the energy decreases
        }
        else if (validity_test == 0){ //Then do regular backtracking
            can_increase = 0;
            alpha *= rho;  // Reduce step size
        }
    }
    return alpha;
}

jcv_real line_search_by_force_local_minima(data* sys, jcv_point* gradient, jcv_point* direction, 
    jcv_real alpha, jcv_real rho, jcv_real c __attribute__((unused))) {
    
    int can_increase = 1;

    jcv_point x_old[sys->N];

    // Store initial positions
    for (int i = 0; i < sys->N; i++) {
        x_old[i] = sys->positions[i];
    }

    // Calculate initial force magnitude squared
    jcv_real initial_force_mag_sq = 0.0;
    for (int i = 0; i < sys->N; i++) {
        initial_force_mag_sq += jcv_lenght_sq(&sys->forces[i]);
    }

    // Calculate initial dot product to ensure we're going downhill
    jcv_real dot_product = 0.0;
    for (int i = 0; i < sys->N; i++) {
        dot_product += jcv_dot(gradient[i], direction[i]);
    }

    // If going uphill, change direction to -gradient (pure gradient descent)
    if (dot_product > 0) {
        dot_product = 0;
        for (int i = 0; i < sys->N; i++) {
            direction[i] = sys->forces[i];
            dot_product += jcv_dot(gradient[i], direction[i]);
        }
    }

    int max_iterations = 20;  // Prevent infinite loop
    int iter = 0;

    while (iter++ < max_iterations) {
        // Take a step
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = x_old[i].x + alpha * direction[i].x;
            sys->positions[i].y = x_old[i].y + alpha * direction[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }

        // Compute new forces
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);

        // Calculate new force magnitude
        jcv_real new_force_mag_sq = 0.0;
        for (int i = 0; i < sys->N; i++) {
            new_force_mag_sq += jcv_lenght_sq(&sys->forces[i]);
        }

        int validity_test = new_force_mag_sq < initial_force_mag_sq;
        if (validity_test && can_increase == 0) {
            break;  // Stop if a Armijo condition is satisfied
        }
        else if (validity_test && can_increase == 1){
            alpha /= rho; //Start with small alpha, increases it while the energy decreases
        }
        else if (validity_test == 0){ //Then do regular backtracking
            can_increase = 0;
            alpha *= rho;  // Reduce step size
        }

    }

    return alpha;  // Return the last tried alpha
}

// Add this function after the existing line search functions but before conjugateGradientStep

jcv_real line_search_forcezero(data* sys, jcv_point* gradient __attribute__((unused)), jcv_point* direction, 
    jcv_real alpha, jcv_real rho __attribute__((unused)), jcv_real c __attribute__((unused))) {
    
    jcv_point x_old[sys->N];
    
    // Constants (similar to LAMMPS)
    const jcv_real ZERO_ENERGY = 1e-12;   // Max allowed energy increase
    const jcv_real GRAD_TOL = 0.1;        // Target reduction for directional derivative
    const jcv_real ALPHA_FACT = 0.1;      // Factor to reduce alpha when backtracking
    const jcv_real MIN_ALPHA = 1e-14;     // Minimum allowed alpha before giving up
    const jcv_real LIMIT_BOOST = 4.0;     // Maximum boost for alpha_del
    
    // Calculate maximum allowed step size (alpha_max)
    jcv_real hmax = 0.0;
    for (int i = 0; i < sys->N; i++) {
        jcv_real h_mag = sqrt(direction[i].x*direction[i].x + direction[i].y*direction[i].y);
        hmax = jcv_max(hmax, h_mag);
    }
    if (hmax < 1e-20) return 0.0;  // If search direction is effectively zero
    
    jcv_real alpha_max = 1.0 / hmax;  // Limit step size based on largest component
    
    // Store initial positions
    for (int i = 0; i < sys->N; i++) {
        x_old[i] = sys->positions[i];
    }
    
    // Calculate initial force projection along search direction (fhCurr)
    jcv_real fhCurr = 0.0;
    for (int i = 0; i < sys->N; i++) {
        fhCurr += jcv_dot(sys->forces[i], direction[i]);
    }
    
    // If fhCurr <= 0, search direction isn't downhill
    if (fhCurr <= 0.0) {
        // Use gradient descent instead
        for (int i = 0; i < sys->N; i++) {
            direction[i] = sys->forces[i];
        }
        fhCurr = 0.0;
        for (int i = 0; i < sys->N; i++) {
            fhCurr += jcv_dot(sys->forces[i], direction[i]);
        }
        if (fhCurr <= 0.0) return 0.0;  // No improvement possible
    }
    
    jcv_real fhOriginal = fhCurr;
    jcv_real engCurr = energy_total(sys);
    jcv_real engOriginal = engCurr;
    
    // Initialize alpha to 0
    alpha = 0.0;
    
    // Initial alpha_del estimation
    jcv_real alpha_init = jcv_max(1e-6, 0.1*fabs(engOriginal)/fhCurr);
    jcv_real alpha_del = jcv_min(alpha_init, 0.5*alpha_max);
    
    // Main line search loop
    int max_iter = 50;  // Prevent infinite loops
    for (int iter = 0; iter < max_iter; iter++) {
        int backtrack = 0;
        jcv_real fhPrev = fhCurr;
        jcv_real engPrev = engCurr;
        
        // Apply increment to alpha
        alpha += alpha_del;
        if (alpha > alpha_max) {
            alpha = alpha_max;  // Could porobably break here
        }
        
        // Take step with new alpha
        sys->N_pbc = sys->N;
        for (int i = 0; i < sys->N; i++) {
            sys->positions[i].x = x_old[i].x + alpha * direction[i].x;
            sys->positions[i].y = x_old[i].y + alpha * direction[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
        
        // Compute new forces and energy
        jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
        sys->sites = jcv_diagram_get_sites(sys->diagram);
        compute_force(sys);
        engCurr = energy_total(sys);
        
        // Compute new directional derivative
        fhCurr = 0.0;
        jcv_real ffCurr = 0.0;  // Not used in this simplified version
        for (int i = 0; i < sys->N; i++) {
            fhCurr += jcv_dot(sys->forces[i], direction[i]);
            ffCurr += jcv_dot(sys->forces[i], sys->forces[i]);
        }
        
        // Check if energy increased
        jcv_real de = engCurr - engPrev;
        if (de > ZERO_ENERGY) {
            backtrack = 1;
        }
        
        // Check if directional derivative is sufficiently reduced
        if (!backtrack && (fabs(fhCurr/fhOriginal) <= GRAD_TOL)) {
            return alpha;  // Success: reduced gradient enough without increasing energy
        }
        
        // Check if we overshot the minimum (derivative changed sign)
        if (fhCurr < 0.0) {
            backtrack = 1;
        }
        
        // Handle backtracking
        if (backtrack) {
            // Move back to previous position
            alpha -= alpha_del;
            for (int i = 0; i < sys->N; i++) {
                sys->positions[i].x = x_old[i].x + alpha * direction[i].x;
                sys->positions[i].y = x_old[i].y + alpha * direction[i].y;
                pbc(&sys->positions[i], sys->L, sys->gamma);
                addBoundary(sys, i);
            }
            
            // Recompute forces at reverted position
            jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
            sys->sites = jcv_diagram_get_sites(sys->diagram);
            compute_force(sys);
            engCurr = engPrev;
            fhCurr = fhPrev;
            
            // Determine new alpha_del
            if (fhCurr < 0.0) {
                // Force changed sign: linearize and solve for zero
                alpha_del *= fhPrev/(fhPrev - fhCurr);
            } else {
                // Energy increased but force didn't change sign
                alpha_del *= ALPHA_FACT;
            }
            
            // Check if alpha_del is too small
            if (alpha_del < MIN_ALPHA) {
                // Restore original positions as fallback
                sys->N_pbc = sys->N;
                for (int i = 0; i < sys->N; i++) {
                    sys->positions[i] = x_old[i];
                    pbc(&sys->positions[i], sys->L, sys->gamma);
                    addBoundary(sys, i);
                }
                jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
                sys->sites = jcv_diagram_get_sites(sys->diagram);
                compute_force(sys);
                return 0.0;  // Failed line search
            }
        } else {
            // Successful step, try to be more aggressive with linearization
            jcv_real boostFactor = LIMIT_BOOST;
            
            // Avoid problems near inflection point
            if (fhPrev > fhCurr) {
                boostFactor = fhCurr/(fhPrev - fhCurr);
            }
            
            boostFactor = jcv_min(boostFactor, LIMIT_BOOST);
            alpha_del *= boostFactor;
            
            // If we've reached alpha_max or directional derivative is close to zero
            if (alpha >= alpha_max || fabs(fhCurr) < 1e-10) {
                return alpha;  // Successful termination
            }
        }
    }
    
    return alpha;  // Return current alpha value after max iterations
}

void conjugateGradientStep(data* sys) {

    //pointer to line_search function  
    jcv_real (*current_line_search)(data*, jcv_point*, jcv_point*, jcv_real, jcv_real, jcv_real);
    #if FORCE_BASED_MINIMIZATION
        current_line_search = line_search_forcezero;
    #else
        #if STAY_IN_LOCAL_MINIMA
        current_line_search = line_search_local_minima;
        #else
        current_line_search = line_search;
        #endif
        int count = 0;
        int minimize_with_force = 0;
        jcv_real old_gnorm = 0.0;
    #endif

    jcv_real new_gnorm = 0.0;

    jcv_diagram_generate(sys->N_pbc, sys->positions, NULL, NULL, sys->diagram);
    sys->sites = jcv_diagram_get_sites(sys->diagram);

    
    compute_force(sys);


    jcv_real chi = 0.0;
    jcv_point gradient[sys->N];
    jcv_point delta[sys->N];

    for (int i = 0; i < sys->N; i++) {
        gradient[i].x = -sys->forces[i].x; // grad(U)
        gradient[i].y = -sys->forces[i].y;
        delta[i] = sys->forces[i]; // -grad(U)
    }

    

    jcv_real alpha = 1;
    

    while (1){

        #if FORCE_BASED_MINIMIZATION
            current_line_search(sys, gradient, delta, fmax(alpha, 1e-3), 0.5, 0.0001);
        #else
            #if STAY_IN_LOCAL_MINIMA
            alpha = current_line_search(sys, gradient, delta, fmax(alpha, 1e-3), 0.5, 0.0001);
            #else
            alpha = current_line_search(sys, gradient, delta, 0.2, 0.5, 0.0001);
            #endif
        #endif


        compute_force(sys);
        new_gnorm = 0;
        for (int i = 0; i < 2 * sys->N; i++) {
            new_gnorm += jcv_lenght_sq(&sys->forces[i]);
        }
        new_gnorm = JCV_SQRT(new_gnorm/sys->N);
        //printf("gnorm = %e, E = %.10e, alpha = %.10e\n", new_gnorm, energy_total(sys), alpha);

        #if FORCE_BASED_MINIMIZATION == 0
            // counter before using force based minimization criteria
            if (jcv_abs(new_gnorm - old_gnorm) < 1e-14){
                count++;
            }
            else{
                count = 0;
            }

            if (count > 3){
                break;
                #if STAY_IN_LOCAL_MINIMA
                current_line_search = line_search_by_force_local_minima;
                #else
                current_line_search = line_search_by_force;
                #endif
                minimize_with_force = 1;
            }
            if (minimize_with_force == 1 && (new_gnorm - old_gnorm) > 0){
                break;
            }
        #endif

        if (new_gnorm < 1e-10){
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
        #if STAY_IN_LOCAL_MINIMA && FORCE_BASED_MINIMIZATION == 0
        old_gnorm = new_gnorm;
        #endif
        //saveTXT(sys);
    }
    
    loggers(sys);
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
            printf("fnorm = %.14lf, vnorm = %lf, P = %lf, dt = %lf, E = %lf\n", fnorm, vnorm, P, dt_now, energy_total(sys));
            if (fnorm > 0.0001){
                alpha_now = 1;
                dt_now = 0.001;
                P = 1;
                vnorm = 1;
                fnorm = 1; 
            }
            
            //saveTXT(sys);
            
            if (count > 200000){
                printf("Failed to converge\n");
                unexpectedClosure(sys);
                exit(3);
            }
            //check_force(sys);
            //exit(3); 
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
        for (int i = 0; i < sys->N; i++){
            sys->positions[i].x += dt_now*sys->velocities[i].x;
            sys->positions[i].y += dt_now*sys->velocities[i].y;
            pbc(&sys->positions[i], sys->L, sys->gamma);
            addBoundary(sys, i);
        }
  

        
    } while (fnorm/sys->L > 1e-8); // L = sqrt((float)N)
    
    loggers(sys);


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

