%% Generate Irregular FEM Mesh
addpath('Meshing files')

%{
This script uses the mesh generation algorithms of Per-Olof Persson and 
Gilbert Strang to generate a triangulation of an irregular domain.

Sources:
-Algorithms: Per-Olof Persson and Gilbert Strang
-Documentation: "A simple mesh generator in MATLAB"
-Available from: SIAM Review, Vol. 46, No. 2, pp. 329?345


Author of this script: Matthew Geleta
Date: 19/03/2017
%}

% Parameters:
h0 = 1/20; % cell diameter
% coordinates and radii of holes:
p1x = -0.45; p1y = 0.45; R1 = 0.25;
p2x = -0.45; p2y = -0.45; R2 = 0.15;
p3x = 0.4; p3y = 0; R3 = 0.4;
Rout = 1; % radius of outer boundary

% Distance function for three circular holes in circular plate:
fd = @(p) ddiff( ddiff( ddiff(dcircle(p,0,0,Rout),dcircle(p,p1x,p1y,R1))...
    ,dcircle(p,p2x,p2y,R2)), dcircle(p,p3x,p3y,R3));

% Mesh generation; coordinates p and triangulation t:
[p,t] = distmesh2d(fd,@huniform,h0,[-1,-1;1,1],[]);

% Display mesh:
shg;