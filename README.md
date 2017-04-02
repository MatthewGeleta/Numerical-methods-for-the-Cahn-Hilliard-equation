# Numerical-methods-for-the-Cahn-Hilliard-equation

This repository is based on the paper "Numerical Methods for the Cahn-Hilliard Equation", by Matthew Geleta, submitted as part of an MSc degree in Mathematics at The University of Oxford.

--------------------------
Run the script "Example_Script.m" for a demonstration of the finite difference solver, and a pretty simulation.

Run the script "FEM_Cahn_Hilliard_Irregular.m" to view a finite element solution to the stead-state
Cahn-Hilliard equation on an irregular domain.

Run the script "FEM_Cahn_Hilliard_Rectangular.m" to view a finite element solution to the stead-state
Cahn-Hilliard equation on a rectangular domain.
--------------------------

This repository includes MATLAB code for:
-Five finite-difference schemes for the two-dimensional Cahn-Hilliard equation with Neumann boundary conditions.
-A finite element scheme for the steady-state Cahn-Hilliard equation.
-Functions to compute area integrals using a two-dimensional Simpson's rule.
-Functions to generate animations showing the evolution of the Cahn-Hilliard system.
-An example script (called "Example_Script") to demonstrate how the finite-difference code works.
-Some background dependency files called internally within some functions.
