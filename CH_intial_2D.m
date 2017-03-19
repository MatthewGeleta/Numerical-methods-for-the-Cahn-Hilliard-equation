function cvecs = CH_intial_2D(N,T,k,seed)
%{
CH_initial_2D initialised the solution matrix cvecs with user specified
initial condition.

Inputs:
-N:     number of grid spacings in one-dimension
-T:     maximum number of time-steps
-k:     specifies which initial condition
-seed:  seeds the random number generator

Output:
-cvecs:     initialised concentration matrix in the format used by CH2D11,
        CH2D12, CH2D13, CH2D14, CH2D15.

Author: Matthew Geleta
Data: 19/03/2017
%}
    
% Intialise concentration matrix
cvecs = zeros((N+1)^2,T);
% Select initial condition
    switch(k)
        case 1 % Random initial condition
            rng(seed);
            cvecs(:,1) = 2*rand((N+1)^2, 1) - 1;
        case 2 % Smooth cosine initial condition
            x = linspace(0,1,N+1); % get mesh spacings
            [X,Y] = meshgrid(x,x); % put on meshgrid
            cmat = cos(2*pi*X).*cos(pi*Y); % initial condition
            cvecs(:,1) = cmat(:); % reshape into vector format (lexicographical ordering)
    end  
end