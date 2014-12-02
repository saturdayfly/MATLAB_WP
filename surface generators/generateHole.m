function [X, Y, Z] = generateHole(max_slope, num_point)
% Generate a surface made of a single trough
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

% Define period
ll = 1;

% Define amplitude as a function of maximum slope Psi
a0 = tand(max_slope)*ll/2/pi;

% Define the grid and compute Z data
[X, Y]  = meshgrid(linspace(0,ll/2,num_point));
R       = sqrt(X.^2+Y.^2);
Z       = a0.*(1-cos(2*pi*R/ll));

Z(R>=ll/2) = NaN;


% Below are the analytic expression to be compared to the Wenzel prewetting
% model outputs

%{

function is : 

z = a0.*(1-cos(2*pi*r/ll))   --> r = (ll/2/pi)*acos(1-z/a0)
with r = sqrt(x²+y²)

maximum slope of surface is maximum of

z' = (2*pi*a0/ll)*sin(2*pi*r/ll)
so tan(psi) = 2*pi*a0/ll   --> this gives the definition of a0 above

SLOPE : 

derivative is : z' = (2*pi*a0/ll)*sin(2*pi*r/ll)    substituting r it gives
                z' = (2*pi*a0/ll)*sin(acos(1-z/a0))

so we have 
    - sigma1th = (2*pi*a0/ll)*sin(acos(1-z/a0))               CHECKED
    - sigma2th = sigma1th.^2                                    CHECKED

CURVATURE :

 - from http://math.stackexchange.com/questions/153371/how-do-i-compute-mean-curvature-in-cylindrical-coordinates

r       = (ll/2/pi)*acos(1-z/a0);
dfdr    = (2*pi*a0/ll)*sin(2*pi*r/ll);
dfdt    = 0;
dfdrdt  = 0;
d2fdr2  = a0*(2*pi/ll)^2*cos(2*pi*r/ll)

num     = dfdr.*(r.^2.*(dfdr.^2+1)) + r.*r.^2.*d2fdr2;
den     = 2*(r.^2.*(dfdr.^2+1)).^(3/2);

meanK   = num./den;                                         CHECKED

ISOLINES :

length of isoline is :

isoline = 2*pi*r;
area    = pi*(ll/2)^2;
Lth = isoline/area; % normalized by area



%}