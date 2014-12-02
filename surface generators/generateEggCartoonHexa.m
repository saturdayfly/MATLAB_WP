function [X, Y, Z] = generateEggCartoonHexa(max_slope,num_point)
% Generate a 3D egg cartoon pattern with hexagonal coordinates
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
% Notes:
%   - THis is not the final versions, Ciro proposed some modifications to
%   get a better symmetry
%   - max_slope is not used so far, we directly define the amplitude a0
%       (to directly compare with Segun data)
%
%
% Here is some background about how we convert the hexagonal coordinates
% into cartesian coordinates
%
%{

Transform between hexagonal and cartesian coordinates

1/ Relationships between unit vectors

We assume (u,v) the unit vector of the hexagonal lattice forming an
angle of 60°, and (x,y) the unit vector of the cartesian lattice.
Relationships between (u,v) and (x,y)  are :

u = x
v = (1/2)x + sqrt(3)/2*y

or

x = u
y = (2/sqrt(3))v - (1/sqrt(3))u

2/ Relationship between coordinates of a point in space

Assuming P a point in the plane with coordinates P(X,Y) = Xx+Yy or
P(U,V)=Uu+Vv, relationships between (U,V) and (X,Y) are :

X = U + (1/2)V
Y = (sqrt(3)/2)V

U = X - (1/sqrt(3))Y
V = (2/sqrt(3))Y

We use the above relationships to define the hexagonal egg cartoon on a
cartesian mesh

%}   

% Define amplitude 
a0 = 4;

% Define periods
Lx = 50;
Ly = 50;

% size of the periodic box
lx = Lx;
ly = 2*(sqrt(3)/2)*Ly;

% Define the cartesian  grid and compute Z data
[X, Y] = meshgrid(linspace(0,lx,num_point),linspace(0,ly,num_point));

% Define the hexagonal coordinates
U = X - (1/sqrt(3)).*Y;
V = (2/sqrt(3)).*Y;

Z = a0-a0.* ...
    cos(2*pi*U/Lx).*...
    cos(2*pi*V/Ly);
                            
