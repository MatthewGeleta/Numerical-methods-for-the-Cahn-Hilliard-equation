function fighand = CH2D_Plot_Evolution(cvecs)
%{
CH2D_Plot_Evolution plots the evolution in time of the concentration
distribution in a two-dimensional Cahn-Hilliard simulations.

Input:
-cvecs: concentration vectors in the format of the output of CH2D11,
CH2D12, CH2D13, CH2D14, and CH2D15.

Output:
-fighand: figure handle to the concentration plot of the evolution.

Author: Matthew Geleta
Date: 19/03/2017
%}

% Obtain parameters from cvecs input
T = size(cvecs,2);
Nsq = size(cvecs,1);
N = sqrt(Nsq) - 1;

fighand = figure();
    for n = 1:T
    Cmat = reshape(cvecs(:,n), [N+1,N+1]); % reshape vector to matrix form
    %contourf(Cmat); % contour plot
    %colormap('gray'); % colour scheme
    pcolor(Cmat);
    shading interp; % interpolation scheme
    %surf(Cmat); % surface plot
    set(gca,'visible','off') % remove axes ticks
    pause(0.01); % change rate of plotting
    end
end