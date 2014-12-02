function [X, Y, Z] = generateRough(num_modes,num_point,power)
% Generate a rough surface
%
% Args:
%   - num_modes: integer, the  number of modes the surface is made of
%   - num_point: integer, number of points for the grid in x and y
%                directions. Width/height of the pattern is 1.
%                So resolution is 1/num_point
%   - power: power law setting the decay of the amplitudes
%   - phase(optional): the random phases. if not specified random phases
%                      are generated.
%
% Returns:
%   - X,Y, Z: N by M matrices with x, y and z values of the surface profile
% 

% Message
fprintf('Generating a rough surface with the following parameters :\n');
fprintf('-num_modes: %d\n-num_points: %d\n-power: %d\n',num_modes,num_point,power);

%{
if any(phase)  % phase data not empty
    fprintf('Random phases passed by user.\n');
else
    fprintf('Random phases generated.\n');
end
%}

% define q vector and amplitudes (power law decay set by 'power')
    [kx, ky] = meshgrid(-num_modes:num_modes);        % modes grid
    amp = 2/num_modes/num_modes;                      % amplitude of first mode
    kmag = sqrt(kx.^2+ky.^2);
    a = amp./kmag.^power;                              % amplitudes

% Set amplitude of mode 0 to 0 (i.e. zero mean value)
    a(num_modes+1,num_modes+1) = 0;
    
% set amplitude to zero if kx OR ky is zero (cf. ciro function)
    a(kx==0) = 0;
    a(ky==0) = 0;

% fix a threshold for the ampmlitude
% *** This may change in future versions
    spfrac = 0.25;
    cutoff = amp/(2*(num_modes*spfrac)^2);
    a(a>=cutoff) = cutoff;

% Phase data
    %{
    if any(phase)  % phase data given as input
        fh = reshape(phase,2*num_modes+1,2*num_modes+1); % load phase
    else
        fh = rand(2*num_modes+1); % generate a random phase 0 < fh < 1
    end
    %}
    fh = rand(2*num_modes+1); % generate a random phase

% Size of periodic box
    lx = 1;
    ly = 1;

% define height function (sum of all modes)
    height = @(x,y) sum(sum(a.*sin( 2*pi*(	(kx.*x./lx) + (ky.*y./ly) + fh))));

% Define the grid and compute Z data
    [X, Y] = meshgrid(linspace(0,1,num_point));
    Z = reshape(cell2mat(arrayfun(height, X(:), Y(:), 'UniformOutput', 0)),...
                num_point,num_point);



