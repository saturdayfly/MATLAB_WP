% Non-Gaussian Surface Generator
%
% The purpose of this script is to generate non-Gaussian random surfaces
% with specified standard deviation, autocorrelation length, skewness and
% kurtosis.
%
% References : 
%
% [1] J.-J. Wu, Simulation of non-Gaussian Surfaces with FFT, Tribology
%     International, 2004, 37, 339-346
% [2] Y. Z. Hu and K. Tonder, Simulation of 3-D Random rough surfaces by 2-D
%     digital filter and Fourier Analysis, Int. J. Mach. Tools Manufact.,
%     1992, 32(1), 83-90
% [3] A. Urzica and S. Cretu, Simulation of the non-gaussian roughness with
%     specified values for the high order moments.

%% TO DO

% In part 1, distinguish spatial coordinates k,l from frequency coordinates
% omegax, omegay
%
%


%% Setup

addpath('C:\Users\Renaud\Desktop\Wenzel Prewetting\_MATLAB_WP\libraries\customplots');

%% PART 1 - Generation of isotropic Gaussian surfaces (Ref 2 part 3.1)
% rough surface generated from low-pass filtering of a random series

    m     = 32;
    n     = 32;
    M     = 64;
    N     = 64;
    Fs    = 64;     % Sampling frequency
    sigma = 1;      % rms roughness
    alpha = 0.2;    % ratio between cut-off and sampling frequency

% input sequence of normal iid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rng(2787483)
    eta = normrnd(0,1,N,M);

% define 2-D ideal low-pass filter with Lanczos window %%%%%%%%%%%%%%%%%%%%

    % Define frequency grid
    % filter size is m by n (32*32, half the sample size)
    k = -m/2+1:1:m/2-1;
    l = -n/2+1:1:n/2-1;
    [U, V] = meshgrid(k,l);
    rho = sqrt(U.^2 + V.^2);

    % Compute filter coefficients
    %(from analytic expression, See annexe for demonstration)
    hLPF_ideal = alpha*besselj(1,2*pi*alpha*rho)./rho; 
    hLPF_ideal(m/2,n/2) = pi*alpha^2;   % amplitude of continuous component

    % Create Lanczos window
    % note : ref [2] uses a parameter m_w = 1.6. I assume it is an exponent
    Lanczos = sinc(rho/15).^1.6;
    Lanczos(rho>15) = 0;

    % Compute fft of filter + window
    HLPF = fftshift(fft2(hLPF_ideal.*Lanczos,31,31));

    % plot the filter + window
    figure; subplot(1,2,1);
    mesh(U,V,hLPF_ideal.*Lanczos); colormap([0 0 0]); grid off;
    %oaxes([0 0 0]);
    xlabel('k'); ylabel('l'); zlabel('h(k,l)');
    subplot(1,2,2);
    mesh(U,V,abs(HLPF)); colormap([0 0 0]); grid off
    xlabel('fx'); ylabel('fy'); zlabel('H(fx,fy)');
    %oaxes([0 0 0]);

% Apply the filter to the input sequence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After the filtering Z does not have unit variance, so it is rescaled

    Z = filter2(hLPF_ideal.*Lanczos,eta);
    Z = Z./std(Z(:)); % scale to get unit variance
    [X, Y] = meshgrid(1:M,1:N);

    % plot
    figure; subplot(1,2,1); surf(X,Y,eta);  % input random sequence
    subplot(1,2,2); surf(X,Y,Z);            % filtered sequence

% plot height distribution

    nbins = 30;
    z_binEdges = linspace(min(Z(:)),max(Z(:)),nbins+1);
    z = ( z_binEdges(1:end-1) + z_binEdges(2:end) ) ./ 2; % bin centers
    p_z = histc(Z(:),[z_binEdges(1:end-1) Inf])'; p_z(end) = [];
    p_z = p_z/trapz(z,p_z); % normalize

    figure; subplot(1,2,1);
    plot(z,p_z,'ob',-4:0.1:4,(2*pi)^(-1/2)*exp(-(-4:0.1:4).^2/2),'-r');
    xlabel('z'); ylabel('p(z)');
    legend('height distribution','normal distribution')
    subplot(1,2,2); normplot(Z(:))
    

%% PART 2 - Generation of Gaussian surface with specified ACF (Ref 2 part 3.2)
% non-isotropic rough surface with specified autocorrelation function

m     = 32;
n     = 32;
M     = 64;
N     = 64;
Fs    = 64;     % Sampling frequency
sigma = 1;      % rms roughness

% Specify Autocorrelation function Rxx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % We Use an exponential ACF
    % Rxx(k,l) = sigma² exp[-2.3[(k/Bx)²+(l/By)²]^(1/2)]
    % With Bx and By the correlation lengths (defined at 10%)
    % Size of the ACF function is m by n ( = PSD and filter size)

    % correlation lengths (non-isotropic)
    Bx = 4;
    By = 12;
    
    Rxx = 
    
    
% Generate 2D input sequence (normal iid)

% Compute power spectral density Sxx from autocorrelation

% Determine the frequency response

% Compute the filter coefficients

% Compute the 2D output sequence


%% ANNEXE

% ANNEXE A - demonstration of hLPF_ideal analytical expression (cf. part 1)
% see http://books.google.pt/books?id=d8FMlHewYp0C&pg=PA82&dq=%22jinc+function%22#v=onepage&q=%22jinc%20function%22&f=false
% see also the Hankel transform on Wikipedia (relation to fourier transform)

%{
The inverse 2D Fourier transform of a circular symmetric function is given
in polar coordinates by : 

h(r) = 2*pi* int_0^inf [  rho * H(rho) * J0(2*pi*rho*r) * drho ]
with J0 the zero order Bessel function, rho the frequency and r spatial distance.

considering a circle function of cutoff frequency rhoc,
its inverse fourier transform is then : 

h(r) = 2*pi* int_0^rhoc [  rho * 1 * J0(2*pi*rho*r) * drho ]

we perform the variable change : psi = 2*pi*rho*r which gives : 

h(r) = 1/(2*pi*r²) * int_0^{2*pi*r*rhoc} [ psi*J0(psi)*dpsi ]

Using the identity : int_0^x [ psi*J0(psi)*dpsi ] = x*J1(x) we finally get :

h(r) = rhoc * J1(2*pi*rhoc*r) / r

THis is roughly what is used in reference [2] (K. Tonder) although
it is not clear how normalized frequency is introduced (note that the above
demonstration consider fourier image and not FFT, so sampling frequency does not
enter into account)
%}

% Apply the filter

z = filter2(myfilter,eta);

