%% Generate simulations with user-specified scheme
%{
This script produced numerical simulations of the Cahn-Hilliard equation
for a specific user-specified scheme.

Author: Matthew Geleta
Date: 19/03/2017
%}

%% Fix parameters for simulations
N = 50; % Number of spatial steps in one-dimension
T = 100; % Number of time-step
mu = 1; % mu = dt/dx^2
seed = 10; % Seed for random number generator

FileLocation = 'Large_Sims'; % Location in which to save files
SchemeNumber = '2'; % Scheme number with which to generate simulations

%% Generate simulations for different epsilon values:
epsilons = [0.005, 0.01, 0.015]; % Surface tension to loop through
eplen = length(epsilons);

funcname = ['CH2D1',SchemeNumber];
fh = str2func(funcname);

for j = 1:eplen
    ep = epsilons(j);
    cvec = fh(N,T,ep,mu,seed,1);
    fname = [FileLocation,'/N',int2str(N),'ep',num2str(ep),'.mat'];
    save(fname, 'cvec','seed', 'N', 'T', 'ep','mu');
end



