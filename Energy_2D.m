function En = Energy_2D(C,ep)
%{
Energy_2D computes the Cahn-Hilliard free energy using a 2-dimensional
version of Simpson's rule for energy computation.

Function takes in concentration matrix C such that C(:,t) is the
lexicographically ordered vector of concentration vector at time t.

Output is an energy vector En such that En(t) is the energy at time t.

Assumes integration domain is [0,1]x[0,1]

Note: C(:,t) must have length  such that sqrt(N) is odd.

%}

% Get dimensions
[N,T] = size(C); % N must be square of an odd number
n = sqrt(N); % Number of points in single column
h = 1/(n-1); % width of single column

% Energy vector
En = zeros(T,1);

% Simpson weights
s = zeros(n,1);
s(1) = 1; s(n) = 1; s(2:2:end-1) = 4; s(3:2:end-2) = 2;
S = s*s'; % Weight matrix
S = reshape(S, [N,1]); % Weight vector
S = S*(h^2/9);

% Differentiation matrices
ee = ones(n,1);
D1 = spdiags([-ee, ee], [-1,1], n, n); % First derivative matrix with central differences
D1(1,2) = 0; D1(n,n-1) = 0; % Neumann BCs
D2D = kron(speye(n), D1) + kron(D1, speye(n)); % 2 dimensional differentiation matrix

% Define the appropriate integrand g(c)= F'(c) = ep^2|grad C|^2
g = @(c) (1 - c.^2).^2/4 + ep^2/2 * (D2D*c).^2/(4*h^2);
    for tt = 1:T
        En(tt) = S'*g(C(:,tt));
    end
    
end