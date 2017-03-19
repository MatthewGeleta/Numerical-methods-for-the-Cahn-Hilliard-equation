%% Generate Rectangular mesh (for testing purposes)
addpath('Meshing files')
%{
This script uses the mesh generation algorithms of Per-Olof Persson and 
Gilbert Strang to generate a triangulation of a rectangular domain.

Sources:
-Algorithms: Per-Olof Persson and Gilbert Strang
-Documentation: "A simple mesh generator in MATLAB"
-Available from: SIAM Review, Vol. 46, No. 2, pp. 329?345


Author of this script: Matthew Geleta
Date: 19/03/2017
%}

% Parameters:
h0 = 1/20; % cell diameter

pfix = [-1,-1;-1,1;1,-1;1,1];

fd = @(p) drectangle(p, -1, 1, -1, 1);

[p,t]=distmesh2d(fd,@huniform,h0,[-1,-1;1,1],pfix);

% Display mesh:
shg;