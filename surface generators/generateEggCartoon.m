function [X, Y, Z] = generateEggCartoon(max_slope,num_point)
% Generate a 3D egg cartoon pattern
%
% Args:
%   - max_slope: numeric, the maximum slope of the surface profile
%   - num_point: integer, number of points for the grid in x and y
%                directions. Width/height of the pattern is 1.
%                So resolution is 1/num_point
%
% Returns:
%   - X,Y, Z: N by M matrices with x, y and z values of the surface profile
% 

% Define amplitude as a function of maximum slope Psi
a0 = tand(max_slope)/(sqrt(2)*2*pi);

% size of the periodic box
lx = 2;
ly = 2;

% Define the grid and compute Z data
[X, Y] = meshgrid(linspace(0,2,num_point));
Z = a0.*(cos(2*pi*X)+cos(2*pi*Y));