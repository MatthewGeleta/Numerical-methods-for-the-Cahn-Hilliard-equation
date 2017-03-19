function J = localJ(coords, U, ep)
% Function localDJ computes the integrals defining DJ over a single finite
% element triangle with vertices located at 'coords', using linear Lagrange
% basis elements.

G = [1,1,1;coords']\[0,0;1,0;0,1]; % Jacobian for integral transformation
T = det([1,1,1;coords'])/2; % Area of triangular cell

J = T*( ep^2 * G*G'... % Integral with derivatives
    - [2,1,1;1,2,1;1,1,2]/12 + ... % Linear part of intergal
    [12*U(1)^2+2*(U(2)^2+U(3)^2+U(2)*U(3))+6*U(1)*(U(2)+U(3)),...
    3*(U(1)^2+U(2)^2)+U(3)^2+4*U(1)*U(2)+2*U(3)*(U(1)+U(2)),...
    3*(U(1)^2+U(3)^2)+U(2)^2+4*U(1)*U(3)+2*U(2)*(U(1)+U(3));
    3*(U(1)^2+U(2)^2)+U(3)^2+4*U(1)*U(2)+2*U(3)*(U(1)+U(2)),...
    12*U(2)^2+2*(U(1)^2+U(3)^2+U(1)*U(3))+6*U(2)*(U(1)+U(3)),...
    3*(U(2)^2+U(3)^2)+U(1)^2+4*U(2)*U(3)+2*U(1)*(U(2)+U(3));
    3*(U(1)^2+U(3)^2)+U(2)^2+4*U(1)*U(3)+2*U(2)*(U(1)+U(3)),...
    3*(U(2)^2+U(3)^2)+U(1)^2+4*U(2)*U(3)+2*U(1)*(U(2)+U(3)),...
    12*U(3)^2+2*(U(1)^2+U(2)^2+U(1)*U(2))+6*U(3)*(U(1)+U(2))]/60);

end