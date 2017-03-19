%% Time 2D schemes
%{
Script to calculate the solution time for each two-dimensional scheme as a
function of the mesh resolution, and to produce plots of the results.

Author: Matthew Geleta
Date: 19/03/2017
%}

%
%% Initialise
%
tvec = zeros(5,1);
NameStart = 'CH2D1';
T = 10; % Maximum number of time-steps
ep = 0.05; % Surface tension parameter
mu = 1; % mu = dt/dx^2
seed = 10;
k = 1;
%}

%% Solution times vs problem size:
%
Nvals = [5,10,20, 40, 80, 100]; % Mesh resolutions to loop over
Nlen = length(Nvals);
tmat = zeros(5,Nlen);

for n = 1:Nlen
    N = Nvals(n);
    for j = 1:5
        fname = [NameStart,int2str(j)];
        fh = str2func(fname);
        [~,t] = fh(N,T,ep,mu,seed,k); % Computation of solution times
        tvec(j) = t;
    end
    tmat(:,n) = tvec;
end
    %}

%% Produce plots of time vs mesh resolution
%
fg1 = figure(1);
%set(fg1, 'Position', [5 5 1000 1000]);
for j = 1:5
    loglog(Nvals,tmat(j,:))
    hold on;
end
    xlim([5,100])
    title('Time scaling with problem size')
    ylabel('Solution time (sec)')
    xlabel('Problem size N')
    legend('a','b','c','d','e')
%}