# Voronoi

This code implements a molecular dynamics simulation of points based on their Voronoi diagram. The code includes functions to generate Voronoi diagrams using [jcv_voronoi](https://github.com/JCash/voronoi), compute forces, and calculate energy based on the Voronoi cells. It can perform energy minimization through FIRE or Conjugate Gradient or overdamped dynamics through RK4, RK45 with an adaptative timestep or euler.

## Compilation

The project uses a Makefile to manage the build process. There are three build configurations available:

1. **Default Configuration**: Optimized for performance.
2. **Debug Configuration**: Includes debugging information.
3. **Full Debug Configuration**: Includes debugging information and sanitizers for detecting memory errors.

### Build Commands

To compile the project, navigate to the project directory and use the following commands:

- **Default Configuration**:
    ```sh
    make

- **Debug Configuration**:
    ```sh
    make debug

- **Full Debug Configuration (with sanitizer)**:
    ```sh
    make full_debug

## Usage

After compiling the project, you can run the generated binary:
```sh
./voronoi
```

This will execute the Voronoi diagram generation and related computations as defined in the voronoi.c file. You can add additional options such as `-N 1000` to change the number of particles (here 1000).

## Dependencies

Ensure you have the following dependencies installed:
- GCC
- Make
- Openmp (optional)
