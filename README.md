# Voronoi

This codes implement a molecular dynamics simulation of points based on their Voronoi diagram. The code includes functions to generate Voronoi diagrams using jcv_voronoi, compute forces, and calculate energy based on the Voronoi cells.

## Files

- `voronoi.c` and `voronoi.h`: Main file that initializes the system and runs the simulation.
- `helper.c` and `helper.h`: Helper functions used throughout the project.
- `force.c` and `force.h`: Functions to compute forces based on the Voronoi cells.

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