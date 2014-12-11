% Tests of ACF, FFT and PSD computation for a 1D signal

close all
clc
clear

% Define spatial signal
dx   = 1e-3;
Fsx  = 1/dx; % sampling frequency
X    = 0:dx:2;
L    = length(X);
Lfft = 2^nextpow2(L);
dfx  = 1/(Lfft*dx); % = Fsx/L   L should be a power of 2 for FFT alrotiyhm to be consistent
% 100 Hz sine signal
T   = 0.1;  % period of the signal
Z   = sin(2*pi*X/T);

% Compute PSD from FFT
fourier_image  = fftshift(fft(X,Lfft)/Lfft);
PSDfromFFT = abs(fourier_image).^2;
freq           = -Fsx/2:dfx:Fsx/2-dfx; % Frequency axis

% Autocorrelation of sine signal
Rxx = xcorr(X,'biased');
lag = (-(L-1):1:(L-1))*dx;

% PSD based on FFT of autocorrelation function
PSDfromRxx = abs(fftshift(fft(Rxx,Lfft)/Lfft));

% Parseval sums
fprintf('\nParseval sums : \n\n');
fprintf('Spatial sum = %f\n', sum(X(:).^2));
fprintf('Frequency sum (from FFT) = %f\n', sum(PSDfromFFT(:))*Lfft);
fprintf('Frequency sum (from ACF) = %f\n', sum(PSDfromRxx(:))*Lfft);

% Plots
figure('units','normalized','outerposition',[0 0 1 1], ...
        'name', 'Signal and Autocorrelation function');
subplot(2,1,1); plot(X,Z,'-b'); title('Signal z(x)');
subplot(2,1,2); plot(lag,Rxx,'-r'); title('Autocorrelation');

figure('units','normalized','outerposition',[0 0 1 1], ...
        'name', 'Power spectral densities');
plot(freq,PSDfromFFT,'-or'); title('PSD based on FFT'); hold on;
plot(freq, PSDfromRxx, '-ob'); title('PSD based on autocorrelation');
legend('PSD based on FFT','PSD based on autocorrelation'); hold off;
set(gca,'yscale','log');