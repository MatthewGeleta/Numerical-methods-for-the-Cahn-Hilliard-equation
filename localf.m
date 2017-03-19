function f = localf(coords, U, ep)
% Function localJ computes the energy functional J over a single finite
% element triangle with vertices located at 'coords', using linear Lagrange
% basis elements.

G = [1,1,1;coords']\[0,0;1,0;0,1]; % Jacobian for integral transformation
T = det([1,1,1;coords'])/2; % Area of triangular cell

f = T* ( (ep^2 * G*G'... % Integral with derivatives
    - [2,1,1;1,2,1;1,1,2]/12)*U + ... % Linear part of intergal
    [4*U(1)^3+ U(2)^3+U(3)^3+3*U(1)^2*(U(2)+U(3))+2*U(1) ... % Nonlinear integral
    *(U(2)^2+U(3)^2)+U(2)*U(3)*(U(2)+U(3))+2*U(1)*U(2)*U(3);
    4*U(2)^3+ U(1)^3+U(3)^3+3*U(2)^2*(U(1)+U(3))+2*U(2) ...
    *(U(1)^2+U(3)^2)+U(1)*U(3)*(U(1)+U(3))+2*U(1)*U(2)*U(3);
    4*U(3)^3+ U(2)^3+U(1)^3+3*U(3)^2*(U(2)+U(1))+2*U(3) ...
    *(U(2)^2+U(1)^2)+U(2)*U(1)*(U(2)+U(1))+2*U(1)*U(2)*U(3)]/60 );

end