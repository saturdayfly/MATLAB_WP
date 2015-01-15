function [X, Y, Z, attributes] = generateEggCartoonHexa(amplitude,phi,num_point)
% Generate a 3D egg cartoon pattern with hexagonal coordinates
% It is possible to vary the saddle point with the phi parameter
%
% Args:
%   - amplitude: numeric, the amplitude (maximum height) of the surface
%   - phase parameter in degree to shift the saddle point up and down
%     (introduced by Ciro)
%   - num_point: integer, number of points for the grid in x and y
%                directions. Dimension of the pattern is 1 by 1 unit.
%                So resolution is 1/num_point
%
% Returns:
%   - X,Y, Z: N by M matrices with x, y and z values of the surface profile
%   - attributes : a structure containing the following variables :
%       .unit: the unit of X,Y and Z coordinates
%       .width: width of the surface in unit (x dimension)
%       .height: height of the surface in unit (y dimension)
%       .x_resolution: resolution in unit per pixel
%       .y_resolution: resolution in unit per pixel
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some background about how we convert the hexagonal coordinates %
% into cartesian coordinates                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%
%% Code %
%%%%%%%%%

% Define amplitude %%%%
    a0 = amplitude;

% Define periods %%%%
    ll = 1;
    Lx = ll;
    Ly = ll;

%                       Previous version commented below
%{
% size of the periodic box %%%%
    lx = Lx;
    ly = 2*(sqrt(3)/2)*Ly;

% Define the cartesian  grid and compute Z data
    [X, Y] = meshgrid(linspace(0,lx,num_point), ...
                      linspace(0,ly,num_point));

% Define the hexagonal coordinates
    U = X - (1/sqrt(3)).*Y;
    V = (2/sqrt(3)).*Y;

    Z = a0-a0.* ...
        cos(2*pi*U/Lx).*...
        cos(2*pi*V/Ly);
%}

%                            Ciro's version below

% size of the periodic box %%%%
    lx = Lx;
    ly = Ly/sqrt(3);

% Define the cartesian  grid and compute Z data
    [X, Y] = meshgrid(linspace(0,lx,num_point), ...
                      linspace(0,ly,num_point));

% Define the hexagonal coordinates
    u1 =  X;
    u2 =  X/2 + Y*sqrt(3)/2;
    u3 = -X/2 + Y*sqrt(3)/2;
    
    Z = a0*(sin(4*pi*u1/ll+phi/2/pi) + ...
            sin(4*pi*u2/ll) + ...
            sin(4*pi*u3/ll)); 
    
    attributes.width        = 1;
    attributes.height       = 1;
    attributes.x_resolution = 1/num_point;
    attributes.y_resolution = 1/num_point;
    attributes.unit         = '1';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appendix : Evolver function from Ciro %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
ll      = 1;            % x size of the box (with pbc)
lx      = ll;           % x size of the box (with pbc)
ly      = ll/sqrt(3);	% y size of the box (with pbc)
h0      = hp/2;         % this is the film thickness

amp     = 0.03;
%amp    = tan(psi*pi/180)*ly/2/pi;


 function real u1 ( real xx, real yy)
   { local uu;
     uu := xx;
      return uu;
   };

 function real u2 ( real xx, real yy)
   { local uu;
     uu := xx/2+yy*sq3/2;
      return uu;
   };

 function real u3 ( real xx, real yy)
   { local uu;
     uu := -xx/2+yy*sq3/2;
      return uu;
   };


// roughness: height of the surface, centered in h0
 function real hr ( real xx, real yy)
   { local lh;
     lh :=h0+amp*(sin(4*pi*u1(xx,yy)/ll+phi/2/Pi)+sin(4*pi*u2(xx,yy)/ll)+sin(4*pi*u3(xx,yy)/ll) ); 
      return lh;
   };
%}