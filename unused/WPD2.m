function output = WPD2(p_f,z,display)

%%  function output = WPD2
%
%   Compute Wetting phase diagram following the method described in :
%
%   [2] S. Herminghauss, Eur. Phys. J. E., 35(43), 2012
%       "Wetting, Spreading and adsorption on randomly rough surfaces"
%
%   Input : 
%       - p_f       : N*1 column vector containing height distribution
%       - z         : 1*N row vector containing z     values (center of bins used to compute p_f)
%       - display   : set to 'on' to turn on display of graphical resolution (otherwise 'off')
%
%   Output : Structure containing the following fields :
%       - T
%       - rho
%       - R
%       - L
%       - W
%       - thetaf
%       - thetaw
%       - Hf
%       - H_vec
%       - theta_vec
%       - AVG_FILM_THICKNESS
%       - VOLUME
%       - Hc

if strcmp(display,'on')
    showres = 1;
elseif strcmp(display,'off')
    showres = 0;
else
    disp('Invalid input argument ''display'' in function WPD1. Display should be ''on'' or ''off''');
end

% Parameters
        
    span = 0.2; % parameter for smoothing (part 2)
    nvalues_theta       = 50;   % number of theta values used for phase diagram computation
    nvalues_H           = 50;   % number of H     values used for phase diagram computation


%% Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             PART 4.1 : Compute T, ACF(T) and polynomial coefs           %           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute T
    T = NaN(size(z)); T(1)=NaN;
    for k = 2:1:length(z)
        T(k) = sqrt(2)*erfinv(2*trapz(z(1:k),p_f(1:k))-1);
    end

    % fit with cubic spline
    cs = spline(z,T);   

    % Compute radial autocorrelation of T(f(z))
    TZ = ppval(cs,Z);                   % Apply T to Z values
    [R, rho] = radial_ACF(X,Y,TZ,400);  % Compute radial ACF
    R(1) = 0;                           % fix R(r=0) to 1  *           
                                        % * Sometimes NaN is returned due to bining ...
    
    % Gaussian fit (to estimate correlation length)
    idx1 = ~isnan(rho) & ~isnan(R);
    f = fit(rho(idx1)',R(idx1)','exp(-(x/T)^2)');
    coefs1 = coeffvalues(f); cl = coefs1(1);

    % polynomial fit (to get Dx coefficients)
    idx2 = ~isnan(rho) & ~isnan(R) & rho<cl;
    polystring = '1 - D2*x^2 - D4*x^4';
    poly = fit(rho(idx2)',R(idx2)',polystring);
    coefs2 = coeffvalues(poly); D2 = coefs2(1); D4 = coefs2(2);

%% Part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 PART 2 : Compute q(T), L, W                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L(h) : length of contour line a heigfht h     (ref [2] equations 36)
% W(h) : wetted sample area                     (ref [2] equations 44)

% FOLLOWING IS A TEST TO compare q(T(f)) and p(f)/T'
% (check it gaves the same result and correspond to Gaussian distribution)
%{
    % define binning

        psi_binEdges = linspace(min(TZ(:)),max(TZ(:)),z_nbins+1);       % bins
        psi = ( psi_binEdges(1:end-1) + psi_binEdges(2:end) ) ./ 2;     % bins center

    % compute histogram + normalize

        q_T = histc(TZ(:),[psi_binEdges(1:end-1) Inf]); q_T(end) = [];
        q_T = q_T/length(TZ(:))/(psi_binEdges(2)-psi_binEdges(1));

    % plot
        figure; plot(psi,q_T,'-bo'); hold on;
        plot(psi,(1/sqrt(2*pi))*exp(-psi.^2/2),'-','color','red','linewidth',2);
        dz = z(2)-z(1);
        plot(ppval(cs,z),p_f./gradient(T,dz)','o','color','green');
        title('q(T) (normalized)'); xlabel('\chi [a.u.]');
        legend( 'q(T) (computed from the height profile)', ...
            '(2\pi)^{-0.5}*exp(-\chi²/2)','p(f)/grad(T)');
%}
% CONCLUSION OF TEST :
% We do have the same result and a gaussian with sigma = 1
% So we can use p_f./gradient(T,dz)' to compute q_T

    LS = sqrt(pi*D2)*q_T;
    WS = 0.5*(1+erf(T/sqrt(2)));
               
%% Part 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 PART 3 : Computation of adsobed liquid                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Coexistence case H = 0 --> thetaf = max(4./gradient(T)*sqrt(D2/pi))
    % (ref [2] eq. 42)
    
        span = 0.2; % parameter for smoothing
        thetaf =  max(smooth(z,4./gradient(T)*sqrt(D2/pi),span,'rloess'));
    
        plot(   z,4./gradient(T)*sqrt(D2/pi),'o', ...
                z,smooth(z,4./gradient(T)*sqrt(D2/pi),span,'rloess'),'-'); 
        
    % Case theta = 0 --> Hf = thetaf * max(L/W)     (ref [1] eq.12) 
    
        Hf = thetaf*max(LS(2:end)./(2*WS(2:end))');
    
        
        % TO BE CONTINUED .....

    %% Part 4
% Assign data to output structure

output.sigma1               = sigma1;
output.sigma2               = sigma2;
output.L                    = L;
output.W                    = W;
output.thetaf               = thetaf;
output.thetaw               = thetaw;
output.Hf                   = Hf;
output.H_vec                = H_vec;
output.theta_vec            = theta_vec;
output.AVG_FILM_THICKNESS   = AVG_FILM_THICKNESS;
output.VOLUME               = VOLUME;
output.Hc                   = Hc;

