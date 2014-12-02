% Test of surface reconstruction using iterative fitting

% Create a surface & display

    % parameters
    num_modes = 1;
    power     = 2;
    num_point = 100;
    
    % define q vector and amplitudes
    [kx, ky] = meshgrid(-num_modes:num_modes); % modes grid
    amp = 2/num_modes/num_modes;               % amplitude of first mode
    kmag = sqrt(kx.^2+ky.^2);
    a = amp./kmag.^power;                      % amplitudes

    % Set amplitude of mode 0 to 0 (i.e. zero mean value)
    a(num_modes+1,num_modes+1) = 0;
    % a(kx==0) = 0;
    % a(ky==0) = 0;

    % generate random phases
    rng(1);
    fh = rand(2*num_modes+1);

    % Size of periodic box
    lx = 1;
    ly = 1;

    % define height function (sum of all modes)
    height = @(x,y) sum(sum(a.*sin( 2*pi*(	(kx.*x./lx) + (ky.*y./ly) + fh))));

    % Define the grid and compute Z data
    [X, Y] = meshgrid(linspace(0,1,num_point));
    Z = reshape(cell2mat(arrayfun(height, X(:), Y(:), 'UniformOutput', 0)),...
                num_point,num_point);

    % Display
    figure; surf(X,Y,Z);

% Iterative fitting
    
    fitted_amplitudes = nan(length(-num_modes:num_modes));

    % Get mode 0
    fitted_amplitudes(num_modes+1,num_modes+1) = mean(Z(:));
    Z = Z - fitted_amplitudes(num_modes+1,num_modes+1);
    figure; surf(X,Y,Z);
    
    % Run curve fitting tool with eq. z = cos(2*pi*x)
    fitted_amplitudes(2,3) = -1.968;
    Z = Z - fitted_amplitudes(2,3)*cos(2*pi*X);
    figure; surf(X,Y,Z);
    
    % Run curve fitting tool with eq. z = cos(2*pi*y)
    fitted_amplitudes(1,2) = 3.698;
    Z = Z - fitted_amplitudes(1,2)*cos(2*pi*Y);
    figure; surf(X,Y,Z);
    
    % Run curve fitting tool with eq. z = cos(2*pi*sqrt(x²+y²))
    fitted_amplitudes(1,3) = -2.842;
    Z = Z - fitted_amplitudes(1,3)*cos(2*pi*sqrt(X.^2+Y.^2));
    figure; surf(X,Y,Z);
    
    
    
    
    