function [X, Y, Z] = generateWave()
% Generate an axisymmetric wave pattern
%
% Returns:
%   - X,Y, Z: N by M matrices with x, y and z values of the surface profile
% 

% Define windows [px] :
width = 736;
height = 480;

% Define resolution [px/m] (use same than from Wyko profilometer): 
dx = 0.815*1e-6;
dy = 0.937*1e-6;

% Define amplitudes
a0 = 10e-6;
a1 = 10e-6;
a2 = 10e-6;

% Define periods
q1 = 100e-6;
q2 = 50e-6;

% Define the grid and Z data
x_vec = -width/2:1:width/2-1;
y_vec = -height/2:1:height/2-1;

[X, Y] = meshgrid(dx*x_vec,dy*y_vec);
R = hypot(X,Y);

Z = a0 + ...
    a1.*cos(2*pi*R/q1) + ...
    a2.*cos(2*pi*R/q2);
