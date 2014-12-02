function slope = computeSlope(X,Y,Z,dt,method)
%  Calculates the slope of the surface Z = f(x,y)
%  
% Args:
% 	- X,Y   : M by N matrices, x and y coordinates (from meshgrid)
%   - Z     : M by N matrix, z coordinate
%   - dt    : size of the moving window (ROI = dt*dt pixels), must be odd
%   - method : character string, specify the method used for fitting
%       - 'fit' : uses the fit function of matlab (with fittype = 'poly11')
%                 this method potentially allow to fit with more complex surface.
%       - 'mldivide' : uses  "\" operator to solve system of nonlinear
%                      equations (fits with a plane). Faster than 'fit'
%       - 'evans' : Use the approach of evans [1] consisting in fitting a
%                   3*3 ROI with a 4th order polynomial. ROI Z indices are:
%
%                   Z1  Z2  Z3
%                   Z4  Z5  Z6
%                   Z7  Z8  Z9
%
%   Returns:
%       - slope : a M by N matrice containing the slope at each point
%
%   References:
%   [1] I. V. FLORINSKY, Accuracy of local topographic variables derived from digital
%       elevation models. Int. J. Geographical Information Science, vol. 12, N. 1 1998
%

% Initialize output with NaN and check if 'dt' is odd

    slope = NaN(size(Z));
    if mod(dt,2)==0 && dt >2
        fprintf('Error in slope_computation, input argument dt must be odd and >2.\n');
        return;
    end

% set parameters relevant to selected computation method

    if strcmp(method,'fit')         % Set up fittype and options.
        ft = fittype( 'poly11' );   % plane fit by default
        opts = fitoptions( ft );
        opts.Lower = [-Inf -Inf -Inf];
        opts.Upper = [Inf Inf Inf];
    elseif strcmp(method,'mldivide')
    elseif strcmp(method,'evans')
        dx = X(1,2)-X(1,1); % resolution
        dy = Y(2,1)-Y(1,1);
    else
        fprintf('Error in slope_computation, unknown method. Choose "fit" or "mldivide".\n')
    end
    
% Compute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
h = waitbar(0,'Initializing waitbar...');

% Calculate min and max values of x,y according to dt
% (we do not consider pixels on the image borders)

    k = (dt-1)/2;
    imin = k+1;
    imax = size(Z,1)-(k+1);
    jmin = k+1;
    jmax = size(Z,2)-(k+1);

for i = imin:1:imax
    waitbar(i/imax,h,sprintf('Computing slopes. %.0f%% along...',i/imax*100))
    for j = jmin:1:jmax
        
        % extract ROI of size dt*dt  (with k = (dt-1)/2)
        ZROI = Z(i-k:i+k,j-k:j+k);
        XROI = X(i-k:i+k,j-k:j+k);
        YROI = Y(i-k:i+k,j-k:j+k);
        
        if strcmp(method,'fit')
            [fitresult, ~] = fit( [XROI(:), YROI(:)], ZROI(:), ft, opts );
            coef = coeffvalues(fitresult);
            ROI_slope = sqrt(coef(2)^2+coef(3)^2);
        elseif strcmp(method,'mldivide')
            Const = ones(size(XROI(:)));            % Vector of ones for constant term
            coef = [XROI(:) YROI(:) Const]\ZROI(:); % Find the coefficients
            ROI_slope = sqrt(coef(1)^2+coef(2)^2);
        elseif strcmp(method,'evans')               % Use evans method cf. ref [1]
            p = (ZROI(1,3) + ZROI(2,3) + ZROI(3,3) - ZROI(1,1) - ZROI(2,1) - ZROI(3,1))/(6*dx);
            q = (ZROI(3,1) + ZROI(3,2) + ZROI(3,3) - ZROI(1,1) - ZROI(1,2) - ZROI(1,3))/(6*dy);
            ROI_slope = sqrt(p^2+q^2);
        end

        slope(i,j) = ROI_slope;
        
    end
end

close(h);