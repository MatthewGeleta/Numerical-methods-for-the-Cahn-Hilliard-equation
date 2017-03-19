%% Finite Element Method:
%   for steady-state Cahn-Hilliard on an irregular domain.

%{
Script implements the finite element method and Newton-Raphson iteration to
find a steady-state solution to the Cahn-Hilliard equation in
two-dimensions.

The foundation of this program is the paper:
"Remarks around 50 lines of Matlab: short finite element implementation",
by Jochen Alberty et al.  (Numerical Algorithms 20.2-3 (1999): 117-137)

Author: Matthew Geleta
Date: 19/03/2017

%}

clear all;

%% Offline computations
%
%% Mesh generation
%
Generate_Irregular_FEM_Mesh;
%}
%% Set parameters
%
MaxIts = 50; % max Newton iterations
ep = 0.1; % epsilon in CH
tol = 10^(-10);  % tolerance for Newton iteration
seed = 10;
N = size(p,1); % number of mesh nodes
M = size(t,1); % number of cells in triangulation
%}

%% Get boundary and interior nodes (circles)
%
AllNodes = 1:N; % List of all nodes

% Nodes for each inner circular boundary:
B1 = (sum( (p-[p1x,p1y]).^2, 2) < R1^2+h0/10 );
B2 = (sum( (p-[p2x,p2y]).^2, 2) < R2^2+h0/10 );
B3 = (sum( (p-[p3x,p3y]).^2, 2) < R3^2+h0/10 );

% Nodes for large circular boundary:
Bout = (sum( p.^2, 2) > Rout^2-h0/10 );

InNodes = (B1|B2|B3); % inner boundary nodes
OutNodes = (Bout); % exterior boundary nodes
BdryVec = (InNodes | OutNodes); % indicator vector fo all boundary nodes
Bdry = AllNodes(BdryVec); % all boundary nodes
IN = AllNodes(~BdryVec); % all interior nodes
%}

%% Visualise boundary nodes
%
figure()
scatter(p(B1,1),p(B1,2))
hold on
scatter(p(B2,1),p(B2,2))
scatter(p(B3,1),p(B3,2))
scatter(p(Bout,1),p(Bout,2))
title('Boundary nodes')
legend('Boundary 1', 'Boundary 2', 'Boundary 3', 'Outer boundary')
%}

%% Boundary-free computation (Neumann conditions)
% Comment this section out to empty the set of boundary nodes
%
Bdry = [];
IN = AllNodes;
%}

%% Coefficient initialisation
%
U = zeros(N,MaxIts); % coefficient vector at each Newton iteration

%U0 = cos(2*pi*p(:,1)).*cos(pi*p(:,2)); % smooth cosine initial condition
%rng(seed);
%U0 = 2*rand(N,1) - 1; % random initial condition
U0 = sign(p(:,2)); % discontinuous initial condition

U(:,1) = U0;

U(Bdry,1) = 0; % Dirichlet boundary conditions
%}

%% Newton iteration
%
for n = 1:MaxIts-1
    % Current Newton iterate
    un = U(:,n);
    % Assemble vectors and matrices
    f = sparse(N,1);
    J = sparse(N,N);
    for j = 1:M
       nodes = t(j,:); % vertices of triangles
       coords = p(nodes,:); % coordinates of vertices
       % Compute J by assembling local contributions
       J(nodes,nodes) = J(nodes,nodes) + localJ(coords, un(nodes), ep);
       % Compute f by assembling local contributions
       f(nodes) = f(nodes) + localf(coords, un(nodes), ep);
    end
    
    % Apply Dirichlet conditions
    W = zeros(N,1);
    W(Bdry) = 0; % Dirichlet
    % Solving one Newton step on iterior nodes
    W(IN) = J(IN,IN)\f(IN);
    % Update Newton iterate
    U(:,n+1) = un - W;
    mbreak = n;
    if norm(W) < tol
        fprintf('\n\nNewton-Raphson converged after %i iterations\n\n', mbreak);
        break
    end
end


%% Plot final Newton iterate
%
n = mbreak;
fg1 = figure(1);
set(fg1, 'Position', [5 5 1000 1000]);
trisurf(t,p(:,1),p(:,2),U(:,n),'facecolor','interp');
%pbaspect([1 1 1])
set(gca,'visible','off');
title('Steady State: FEM with linear Lagrange basis');
zlabel('Concentration');
%}

%% Plot evolution of Newton iterates
% Uncomment this section to view the convergence of the Newton iteration
%{
figure()
for n = 1:mbreak
    trisurf(t,p(:,1),p(:,2),U(:,n),'facecolor','interp')
    %pbaspect([1 1 1])
     %set(gca,'visible','off');
     pause(0.2)
end
%}