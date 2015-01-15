clear all;
close all;

%% Load a surface

[X, Y, Zdata, attributes] = importFromGwyddion();
Z = Zdata{1};

% spatial variables
dx      = attributes.x_resolution; % x resolution
dy      = attributes.y_resolution; % y resolution
[Lx,Ly] = size(Z);                 % Sample size in pixel, x and y directions
Lxfft   = 2^nextpow2(Lx);          % round Lx to the next power of 2 for FFT
Lyfft   = 2^nextpow2(Lx);          % round Ly to the next power of 2 for FFT

% frequency variables
Fsx     = 1/dx;                    % spatial sampling frequency, x direction
Fsy     = 1/dy;                    % spatial sampling frequency, y direction
dfx     = 1/(Lx*dx);               % frequency resolution, x direction
dfy     = 1/(Ly*dy);               % frequency resolution, y direction
    
%% Compute autocorrelation of Z (data goes into structure 'autocorr')
%    - Rxx is the autocorrelation matrix (size 2M-1 by 2N-1)
%    - RxxRadial is the radial average
%    - rho is the vector of interpolated lag values for RxxRadial
%    - XLAG and YLAG are the lag grids (for plotting)

% TO DO : use fft2(Z, number_of_point) !!!!!!!!!!!!
% with number of point a power of 2
% dfx and dfy have to be defined accordingly i.e. Fsx or Fsy divided by
% number of points
%

drho = sqrt(dx^2+dy^2);
[autocorr.Rxx, autocorr.RxxRadial, autocorr.rho] = ...
    computeRadialACF(X,Y,Z,drho,'biased');

% define lag grid
[autocorr.XLAG, autocorr.YLAG] = meshgrid((-(Ly-1):1:(Ly-1))*dx , ...
    (-(Lx-1):1:(Lx-1))*dy);

% Plot
figure; subplot(1,2,1);
plot(rho,R); title('Radial average of surface autocorrelation');
subplot(1,2,2);
surf(autocorr.XLAG, autocorr.YLAG, autocorr.Rxx, 'Edgecolor', 'none');
view(3);

   
%% Compute PSD from autocorrelation

psd.PSDfromRxx = abs(fftshift(fft2(Rxx))); % NORMALISE ?

[psd.xfreq, psd.yfreq] = meshgrid( -Fsx/2:dfx:Fsx/2-dfx , ...
                                   -Fsy/2:dfy:Fsy/2-dfy);
    
%% Compute PSD from surface 2D FFT
% for details regarding normalization of the fft/psd see :
% [1] S. Kandhasamy, "DFT, PSD and MATLAB" (Pdf December 14, 2010,
%     url = http://www.ligo.caltech.edu/~ethrane/Resources/signals/shivaraj_DFT.pdf)
% [2] G. Heinzel, A. Rudiger and R.Schilling, "Spectrum and spectral density
%     estimation by the Discrete Fourier transform(DFT), including a comprehensive
%     list of window functions and some new fat-topwindows".
%     url = http://holometer.fnal.gov/GH_FFT.pdf
% [3] Christopher J. Walsh, Achim J. Leistner, and Bozenko F. Oreb,
%     "Power spectral density analysis of optical substrates for gravitational-wave
%     interferometry", url = http://www.ligo.caltech.edu/~hiro/docs/CSIRO_LIGO_mirror.pdf
% [n] an option is also to scale with dx the fft and with N*df the inverse fft
%     in that case Parseval theorem compare the two integrals sum(z²*dx*dx) and sum(PSD*df*df)
%     see for example http://blaketools.blogspot.de/2012/06/testing-1-2-3.html

% TODO
%
% - go through ref 1
% - finish the code below
% - compute fft of the same size from the autocorrelation
% - compare the results



% Select windows type
w_type = @flattopwin;
W = window(w_type,Lx)*window(w_type,Ly)';

% compute fft
fourier_image = fftshift(fft2(Z));
fourier_image = fourier_image/(Lx*Ly);  % normalize see ref. 1

% compute psd
psd.PSDfromFFT = abs(fourier_image).^2;

% define frequency grid
[psd.xfreq, psd.yfreq] = meshgrid( -Fsx/2:dfx:Fsx/2-dfx , ...
                                   -Fsy/2:dfy:Fsy/2-dfy);

% Check Parseval sums (cf. ref [1])
sum2D = @(a) sum(reshape(a,1,[]));             % sum elements in 2D matrix

spatialSum = sum2D(abs(Z).^2)*dx*dy;           % g units:  [Amplitude/sqrt(m)]
freqSum    = sum2D(psd.PSDfromFFT)*dfx*dfy;    % G units:  [Amplitude*sqrt(m)]

fprintf('Spatial Sum = %1.2f\n', spatialSum);
fprintf('Frequency Sum = %1.2f\n', freqSum);
fprintf('RMS value is %1.2f unit²\n', sqrt(spatialSum/Xres/Yres))

















%% Compare the 2 PSD

figure; subplot(1,2,1);
surf(psd.xfreq, psd.yfreq, psd.PSDfromRxx, 'Edgecolor', 'none');
title('PSD computed from autocorrelation'); view(3);
surf(psd.xfreq, psd.yfreq, psd.PSDfromFFT, 'Edgecolor', 'none');
title('PSD computed from 2D FFT'); view(3);

    
    
    RMS = sqrt(sum((psd.PSDfromFFT(:))));
    sqrt(mean((Z(:).^2)))
    



fourrier_image = fft2(Z{1});
newZ = ifft2(fourrier_image);

figure; surf(X,Y,Z{1});
figure; surf(X,Y,newZ);

fprintf('Properties of input surface :\n')
fprintf('Mean = %f\n', mean(Z{1}(:)))
fprintf('RMS = %f\n', std(Z{1}(:)))
fprintf('skewness = %f\n', skewness(Z{1}(:)))
fprintf('kurtosis = %f\n', kurtosis(Z{1}(:)))
fprintf('Properties of transformed surface :\n')
fprintf('Mean = %f\n', mean(newZ(:)))
fprintf('RMS = %f\n', std(newZ(:)))
fprintf('skewness = %f\n', skewness(newZ(:)))
fprintf('kurtosis = %f\n', kurtosis(newZ(:)))

