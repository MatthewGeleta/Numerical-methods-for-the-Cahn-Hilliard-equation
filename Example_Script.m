%% Example script
%{
This script performs an example computation for the solution and
animation of a two-dimensional Cahn-Hilliard evolution.

Author: Matthew Geleta
Date: 19/03/2017
%}

% Set simulation parameters
N = 60; T = 200; ep = 0.01; mu = 1; seed = 10; k = 1; %random initial condition
% Run one of the schemes
cvecs = CH2D12(N,T,ep,mu,seed,k);
% Plot the evolution
CH2D_Plot_Evolution(cvecs);
% Calculate energy
En = Energy_2D(cvecs,ep);
% Plot energy
figure(); plot(En); title('Energy evolution'); ylabel('Energy'); xlabel('Time-step')
