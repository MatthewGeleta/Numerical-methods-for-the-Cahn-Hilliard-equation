function [cvecs, t] = CH2D13(N,T,ep,mu,seed,k)   
%{
CH2D13 computes the solution two two-dimensional Cahn-Hilliard equation with Neumann
boundary conditions using scheme (c) in "Numerical Methods for the
Cahn-Hilliard Equation", by Matthew Geleta.

Inputs:
-N:     Number of grid-spacings in one-dimension
-T:     Number of time-steps
-ep:    epsilon
-mu:    mu = dt/dx^2
-seed:  seed for random number generator
-k:     case for initial condition

Outputs:
-cvecs: (N+1)^2-by-M matrix of lexicographically ordered concentrations
-t: total time for matrix solutions

%}
    %% Initialise
    %
    h = 1/N;
    %}
    %% Initialial conditions
    %
    [cvecs] = CH_intial_2D(N,T,k,seed);
    %}
    %% Generate constant matrices P, Q, S in A = [P Q; R S]
    %
    [P,Q,S,D] = Generate_2D_Matrices(N,mu);
    %}
    
    %% Solve system
    %
    tic;
    for n = 2:T
        % Load old concentration
        co = cvecs(1:(N+1)^2,n-1);
        % Update nonlinear part
        R = 3*spdiags(co.^2,0, (N+1)^2, (N+1)^2) - (ep^2/h^2)*D - speye((N+1)^2);
        % Form solution matrix
        A = [P,Q; R, S];
        % Update RHS using old concentration
        b = [co; 2*co.^3];
        % x = [cn;wn] is the new iterate
        x = A\b;
        % Get new concentration
        cvecs(:,n) = x(1:(N+1)^2);
    end
    t = toc;
end