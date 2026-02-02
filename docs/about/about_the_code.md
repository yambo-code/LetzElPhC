**LetzElPhC** is a C code designed to compute electron-phonon coupling matrix elements from the outputs of standard Density Functional Theory (DFT) and Density Functional Perturbation Theory (DFPT) calculations.

Currently, it only supports the **Quantum Espresso** code, with long-term plans to support the Abinit code.

The main objective of this project is to facilitate electron-phonon related calculations within the **YAMBO** code (version 5.2 and above), and it works only with **norm-conserving pseudo-potentials**.

The code is released under the [**MIT license**](license.md) and hosted on GitHub: [LetzElPhC GitHub Repository](https://github.com/yambo-code/LetzElPhC).

## Main Features

*   **Symmetry Awareness**: Utilizes full crystal symmetries, ensuring compatibility with the YAMBO code without encountering phase issues.
*   **Parallelization**: Implements multiple levels of parallelization, including:
    *   OpenMP
    *   Plane-wave parallelization
    *   k-point parallelization
    *   q-point parallelization
*   **Efficient I/O**: Utilizes fully parallel I/O via parallel NetCDF-4/HDF5 libraries.
*   **Portability**: Highly portable. The code can be compiled on various CPU architectures and operating systems with minimal to no changes in the source code.
