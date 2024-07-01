## Riemannian and Euclidian augmented Lagrangian algorithms

This directory contains the codes related to the numerical experiments from the article:

R. Andreani, K. R. Couto, O. P. Ferreira, G. Haeser, and L. F. Prudente, Constraint qualifications for a special class of   nonlinear optimization problems and the strong global convergence properties of an augmented Lagrangian method for embedded submanifold, 2024.

- __RiemannianAL__: Contains the codes for the Riemannian augmented Lagrangian algorithm (Riemannian-AL)
- __EuclidianAL__: Contains the codes for the Euclidean augmented Lagrangian algorithm (Euclidean-AL)

Folders RiemannianAL and EuclideanAL include the following routines:

- _main.m_: Main routine to choose the problem to solve
- _evalf.m_: Objective function
- _evalg.m_: Euclidean gradient of the objective function
- _evalcc.m_: Constraints
- _evalnc.m_: Euclidean gradient of the constraints
- _auglag.m_: The augmented Lagrangian algorithm

#### Installing dependencies

Go to the main folder and type:

```sh
install_dependencies
```

#### Instruction

To run the Riemannian-AL algorithm:

```sh
cd RiemannianAL
main
```

To run the Euclidean-AL algorithm:

```sh
cd EuclidianAL
main
```

####  Third-party codes
The directory also contains third-party free codes:

1) __Manopt 7.1__
    - Matlab toolbox for optimization on manifolds
    - Author: Nicolas Boumal
    - https://www.manopt.org/
    - License: GNU General Public License (GPL) version 3
2) __ASA_CG 3.0__
    - Active set algorithm for solving a bound constrained optimization problem
    - Authors: W. W. Hager and H. Zhang
    - https://people.clas.ufl.edu/hager/software/
    - License: GNU General Public License (GPL) version 3
3) __ASA_CG_matlabWrapper 1.3.0.0__
    - Mex interface for bound constrained optimization via ASA
    - Author: Stephen Becker
    - https://www.mathworks.com/matlabcentral/fileexchange/35814-mex-interface-for-bound-constrained-optimization-via-asa
    - License: BSD 2-clause

Please consider the specific licenses of the third-party codes before modifying and/or redistributing them. License information can be found in the code comments or in separate text files in the corresponding directories.

#### License

GNU General Public License (GPL) version 3
