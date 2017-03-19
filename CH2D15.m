function [cvecs, t, itvec] = CH2D15(N,T,ep,mu,seed,k)   
%{
CH2D15 computes the solution two two-dimensional Cahn-Hilliard equation with Neumann
boundary conditions using scheme (e) in "Numerical Methods for the
Cahn-Hilliard Equation", by Matthew Geleta.

Newton's method is used to solve the nonlinear equations that arise at each
step.

Inputs:
-N:     Number of grid-spacings in one-dimension
-T:     Number of time-steps
-ep:    epsilon
-mu:    mu = dt/dx^2
-seed:  seed for random number generator
-k:     case for initial condition

Outputs:
-cvecs: (N+1)^2-by-M matrix of lexicographically ordered concentrations
-t:     total time for matrix solutions
-itvec: vector of interation numbers required at each time-step

%}
    %% Initialise
    %
    h = 1/N;
    MaxIts = 25; % Maximum number of Newton iterations
    tol = 10^(-6); % Tolerance for Newton's method
    itvec = zeros(T,1);
    %}
    %% Initialial conditions
    %
    [cvecs] = CH_intial_2D(N,T,k,seed);
    w = zeros((N+1)^2,T);
    %}
    %% Generate constant matrices P, Q, S in A = [P Q; R S]
    %
    [~,~,~,D] = Generate_2D_Matrices(N,mu);
    %}
    %% Define function for which to compute roots
    %
    % Vector function f = [f1;f2]
    f1 = @(cn,co,wn) cn - co - mu*D*wn;
    f2 = @(cn,co,wn) wn - cn.^3 + cn + (ep^2/h^2)*D*cn;
    f = @(cn,co,wn) [f1(cn,co,wn); f2(cn,co,wn)];

    tic;
        % Compute c(n) from c(n-1) at n via Newton iteration
    for n = 2:T
        % Initial guess
        co = cvecs(:,n-1); cn = co;
        wo = w(:,n-1); wn = wo;
        
        % Commence Newton iteration
        niter = 0; err = tol + 1;
        while(err > tol && niter < MaxIts)
            % Submatrix in Jacobian
            Dn = (ep^2/h^2)*D - spdiags(3*cn.^2,0,(N+1)^2,(N+1)^2)...
                + speye((N+1)^2);
            % Jacobian J
            J = [speye((N+1)^2), -mu*D; Dn, speye((N+1)^2)];
            % Newton step
            s = -J\f(cn,co,wn);
            % Next Newton approximation
            cnn = cn + s(1:(N+1)^2); 
            wnn = wn + s((N+1)^2 + 1:end);
            
            % Calculate error using infinity norm
            err = norm([cnn;wnn]-[cn;wn],'inf');
            % Rename variables
            cn = cnn; wn = wnn;
            niter = niter + 1;
        end
        % Store required number of iterations
        itvec(n) = niter;
        % Update time-step
        cvecs(:,n) = cn; w(:,n) = wn;
    end
    t = toc;
end