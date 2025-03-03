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

- **Debug Configuration**:
    ```sh
    make debug

- **Full Debug Configuration**:
    ```sh
    make full_debug

### Usage

Here is a README.md file for your Voronoi project:

```markdown
# Voronoi

This project implements a Voronoi diagram generator and related computations. The code includes functions to generate Voronoi diagrams, compute forces, and calculate energy based on the Voronoi cells.

## Files

- `voronoi.c`: Main file that initializes the system and runs the simulation.
- `helper.c` and `helper.h`: Helper functions used throughout the project.
- `force.c` and `force.h`: Functions to compute forces based on the Voronoi cells.
- `voronoi.h`: Header file containing definitions and structures used in the project.

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
  ```

- **Debug Configuration**:
  ```sh
  make debug
  ```

- **Full Debug Configuration**:
  ```sh
  make full_debug
  ```

### Clean Up

To remove the compiled binaries, use the following command:
```sh
make clean
```

## Usage

After compiling the project, you can run the generated binary:
```sh
./voronoi
```

This will execute the Voronoi diagram generation and related computations as defined in the voronoi.c file.

## Dependencies

Here is a README.md file for your Voronoi project:

```markdown
# Voronoi

This project implements a Voronoi diagram generator and related computations. The code includes functions to generate Voronoi diagrams, compute forces, and calculate energy based on the Voronoi cells.

## Files

- `voronoi.c`: Main file that initializes the system and runs the simulation.
- `helper.c` and `helper.h`: Helper functions used throughout the project.
- `force.c` and `force.h`: Functions to compute forces based on the Voronoi cells.
- `voronoi.h`: Header file containing definitions and structures used in the project.

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
  ```

- **Debug Configuration**:
  ```sh
  make debug
  ```

- **Full Debug Configuration**:
  ```sh
  make full_debug
  ```

### Clean Up

To remove the compiled binaries, use the following command:
```sh
make clean
```

## Usage

After compiling the project, you can run the generated binary:
```sh
./voronoi
```

This will execute the Voronoi diagram generation and related computations as defined in the voronoi.c file.

## Dependencies

Ensure you have the following dependencies installed:
- GCC (GNU Compiler Collection)
- Make

## License

This project is licensed under the MIT License.
```

This README provides a brief description of the project, instructions on how to compile it using the Makefile, and usage information.
This README provides a brief description of the project, instructions on how to compile it using the Makefile, and usage information.