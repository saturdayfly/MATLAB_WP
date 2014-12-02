function [X, Y, Z] = generateRough2D(num_modes,num_point,rep)
% Generate a rough surface in 2D
%
% Args:
%   - num_modes: integer, the  number of modes the surface is made of
%   - num_point: integer, number of points for the grid in x and y
%                directions. Width/height of the pattern is 1.
%                So resolution is 1/num_point
%   - rep: integer, number of repetition of the pattern
%
% Returns:
%   - X,Y, Z: N by M matrices with x, y and z values of the surface profile
% 

% Random phases given by Segun :
phase_vector = [0.2168889046 3.049983978 2.985641003 -0.8974382877 -2.321329594 ...
3.249774933 -0.4285540581 -2.15332365 -2.540319681 1.710004807 -0.3272564411 ...
0.8007757664 0.6551880836 -2.523090363 2.844948769 3.818153858 -2.557318211 ...
3.937945366 1.374811649 -2.283641577];

% Parameters

% hh      = 1;  %Film thickness, info only
L_x     = 512;  % period
g_rms   = 5;    % rms roughness (initially 5, now 10 following Ciro's input)
g0      = 10;   % vertical shift to get positive values (*1.5 following Ciro's input)
H       = 0.1;  % Hurst exponent


mode        = 1:num_modes;    % mode vector
sum_i       = sum(arrayfun(@(x) x^-(2.0*H+1),1:num_modes));
amplitude   = sqrt((2*g_rms^2/sum_i)*(mode.^-(2.0*H+1)));   % amplitude vector
phase       = phase_vector(1:num_modes); % phase vector

g = @(x,y) g0 + sum(amplitude.*cos(2*pi*mode*x/L_x+phase));

% Define the grid and compute Z data
xnpt = rep*num_point;
ynpt = num_point;
[X, Y]  = meshgrid( linspace(0,rep*L_x,xnpt),linspace(0,L_x,ynpt));
Z = reshape(cell2mat(arrayfun(g, X(:), Y(:), 'UniformOutput', 0)),ynpt,xnpt);
