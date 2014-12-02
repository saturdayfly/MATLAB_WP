clear all;
close all;

% Load a surface
    [X, Y, Z, attributes] = importFromGwyddion();

    dx   = attributes.x_resolution; % x resolution
    dy   = attributes.y_resolution; % y resolution
    Xres = attributes.width;        % size - x direction
    Yres = attributes.height;       % dize - y direction

% Compute radial autocorrelation function (assumes isotropic surface)
% Still have to implement the unbiased normalization
% And check the output with Gwyddion

    drho = sqrt(dx^2+dy^2);
    [R, rho] = radial_ACF(X,Y,Z{1},drho,'biased');
    figure;
    plot(rho,R);
   
% CONTINUE HERE  AND COMPARE FFT(R) TO THE PSD FROM FFT(Z)
% Compute Power spectral density

    df = 1/max(rho);
    f_max = 1/drho;
    freq = 0:df:f_max;
    
    psd = fft(R,length(freq));



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

