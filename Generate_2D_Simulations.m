%% Generate 2D simulations for each scheme
%{
This script generates simulations for the two-dimensional Cahn-Hilliard
equation by looping through the five finite-difference schemes presented.

Author: Matthew Geleta
Date: 19/03/2017
%}


%% Fix parameters for simulations
N = 50; % Number of spatial steps in one-dimension
T = 50; % Number of time-step
ep = 0.05; % Surface tension parameter
mu = 1; % mu = dt/dx^2
seed = 10; % Seed for random number generator

FileLocation = 'C_Arrays'; % Folder name in which to save simulation data

% Cell array for concentrations
CArry = cell(5,1);

%% Generate solutions
for j = 1:5
    funcname = ['CH2D1',int2str(j)];
    fh = str2func(funcname);
    cvec = fh(N,T,ep,mu,seed,1);
    CArray{j} = cvec;
    clear('cvec');
end
%}

%% Save data
%
fname = [FileLocation,'/C_array_2D_ep',num2str(ep),'.mat'];
save(fname, 'CArray','seed', 'N', 'T', 'ep','mu')
%}