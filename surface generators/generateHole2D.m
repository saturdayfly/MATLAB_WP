function [X, Y, Z] = generateHole2D(max_slope,num_point,rep)
% Generate a 3D egg cartoon pattern
%
% Args:
%   - max_slope: numeric, the maximum slope of the surface profile
%   - num_point: integer, number of points for the grid in x and y
%                directions. Width/height of the pattern is 1.
%                So resolution is 1/num_point
%   - rep: number of repetition of the pattern
%
% Returns:
%   - X,Y, Z: N by M matrices with x, y and z values of the surface profile
% 

% Define period
ll = 1;

% Define amplitude as a function of maximum slope Psi
a0 = tand(max_slope)*ll/2/pi;

% Define the grid and compute Z data
[X, Y]  = meshgrid( linspace(-rep*ll/2,rep*ll/2,rep*num_point),linspace(-ll/2,ll/2,num_point));
Z       = a0.*(1-cos(2*pi*X/ll));

% Some analytical background

%{

function is : 

z = a0.*(1-cos(2*pi*x/ll))   --> x = (ll/2/pi)*acos(1-z/a0)

maximum slope of surface is maximum of

z' = (2*pi*a0/ll)*sin(2*pi*x/ll)
so tan(psi) = 2*pi*a0/ll   --> this gives the definition of a0 above

SLOPE : 

derivative is : z' = (2*pi*a0/ll)*sin(2*pi*x/ll)    substituting x it gives
                z' = (2*pi*a0/ll)*sin(acos(1-z/a0))

so we have 
    - sigma1th = (2*pi*a0/ll)*sin(acos(1-z/a0))                 CHECKED
    - sigma2th = sigma1th.^2                                    CHECKED

CURVATURE :

x       = (ll/2/pi)*acos(1-z/a0);
dfdx    = (2*pi*a0/ll)*sin(2*pi*x/ll);

d2fdx2  = a0*(2*pi/ll)^2*cos(2*pi*x/ll)

num     = d2fdx2;
den     = (1+dfdx.^2).^(3/2);

meanK   = num./den/2;       factor 1/2 coz its mean curvature   CHECKED

ISOLINES :

length of isoline is :

isoline = 2;
area    = 1;
Lth = isoline/area; % normalized by area



%}