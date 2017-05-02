function [P,Q,S,D] = Generate_2D_Matrices(N,mu)
%{
Generate_2D_Matrices creates sparse square matrices for linear equation
solvers in the Cahn-Hilliard schemes CH2D11, CH2D12, and CH2D13.

Inputs:
-N: number of mesh spacings in one-dimension
-mu: mu = dt/dx^2

Outputs:
-Matrices to be used in a linear solver.

Authors: Matthew Geleta and Conor McMeel.
Date: 19/03/2017
%}

    P = speye((N+1)^2);
    S = -speye((N+1)^2);
    % One dimensional Laplacian (second derivative matrix) without boundary
    % conditions
    ee = ones((N+1)^2, 1);
    fd1D = spdiags([ee, -2*ee, ee], (-1:1), N+1, N+1);
    % Two dimensional Laplacian
    fdM = kron(speye(N+1), fd1D) + kron(fd1D, speye(N+1));
    D = fdM;
    % Impose Neumann boundary conditions
    for i = 1:(N+1)
       D(i, i+1+N) = 2;
       D((N+1)^2-i+1, (N+1)^2-i-N) = 2;
       D((i-1)*(N+1)+1, (i-1)*(N+1)+2) = 2;
       D((i-1)*(N+1)+N+1, (i-1)*(N+1)+N) = 2; 
    end
    Q = -mu*D;
end