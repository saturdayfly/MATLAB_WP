function [X, Y, Z, attributes] = generateRough3D(num_modes,num_point,varargin)

% q vector and amplitude
[kx, ky] = meshgrid(-num_modes:num_modes); % modes grid
amp = 0.1/2/num_modes;                     % amplitude of first mode
a = amp./sqrt(1+kx.^2+ky.^2);              % amplitudes grid with power law 1/|k|

% Phase data
if ~isempty(varargin)                   %   phase data given as input
    fh = reshape(varargin{1},2*num_modes+1,2*num_modes+1); %   load phase from .dat file
else
    fh = rand(2*num_modes+1);                  % generate a random phase 0 < fh < 1
end

% note :    The mean height of the surface will be : h_mean = amp*sin(fh(21,21)*2*pi)
%           It is fixed by "amp" only if fh(k==0) == 0  (not the case here)

% set amplitude to zero if kx OR ky is zero (cf. ciro function)
a(kx==0) = 0;
a(ky==0) = 0;
a(kx==0 & ky==0) = 0; % this set the mean elevation

% Size of periodic box
lx = 1;
ly = 1;

% height function (sum of all modes)
height = @(x,y) sum(sum(a.*sin( 2*pi*(	(kx.*x./lx) + (ky.*y./ly) + fh))));

% Define the grid and compute Z data
[X, Y] = meshgrid(linspace(0,1,num_point));
Z = reshape(cell2mat(arrayfun(height, X(:), Y(:), 'UniformOutput', 0)),num_point,num_point);


    attributes.width        = 1;
    attributes.height       = 1;
    attributes.x_resolution = 1/num_point;
    attributes.y_resolution = 1/num_point;
    attributes.unit         = '1';