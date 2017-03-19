function Mass = Mass_2D(C)
%{
Mass_2D computes the Cahn-Hilliard mass using a 2-dimensional
version of Simpson's rule for energy computation.

Function takes in concentration matrix C such that C(:,t) is the
lexicographically ordered vector of concentration vector at time t.

Output is an energy vector Mass such that Mass(t) is the mass at time t.

Assumes integration domain is [0,1]x[0,1]

Note: C(:,t) must have length  such that sqrt(N) is odd.

%}

% Get dimensions:
[N,T] = size(C); % N must be square of an odd number
n = sqrt(N); % Number of points in single column
h = 1/(n-1); % width of single column

% Mass vector
Mass = zeros(T,1);

% Simpson weights
s = zeros(n,1);
s(1) = 1; s(n) = 1; s(2:2:end-1) = 4; s(3:2:end-2) = 2;
S = s*s'; % Weight matrix
S = reshape(S, [N,1]); % Weight vector
S = S*(h^2/9);

    for tt = 1:T
        Mass(tt) = S'*C(:,tt);
    end
end