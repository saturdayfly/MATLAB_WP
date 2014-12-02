function outputs = fftAnalysis(X,Y,Z)

%	Script to perform Fourier analysis of a rough surface
%	Written by R. DUFOUR - MPIDS/DMI - renaud.dufour@ds.mpg.de
%	Version : September 2014
%
%   This code is to be optimized (so far no interpolation of the polar coordinates)
%   With something like :
%
%    sampled_radial_slice = interp2(X,Y,imgfp,x,y);
%    radial_average(radius) = mean(sampled_radial_slice);
%
%   INPUTS :
%       - X,Y       : M by N matrices containing x and y coordinates (typically obtainned from meshgrid)
%       - Z         : M by N matrice containing z values
%
%   OUTPUTS :
%       - 
%       - 
%       - 
%
%
% For info 
%    
%    Direct space    |	Frequency space
%    -------------------------------
%    Pixel size dx   |   Image size 2/dx
%    Pixel size dy   |   Image size 2/dy
%    Image size Lx   |   Pixel size 2/Lx
%    Image size Ly   |   Pixel size 2/Ly
%    Resolution Nx   |   Resolution Nx
%    Resolution Mx   |   Resolution Mx
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [N, M] = size(Z);       % Size of input data [px]
    Xres = range(X(:));     % Real size [unit]
    Yres = range(Y(:));     % Real size [unit]
    
    Fs_x = M/Xres;           % Pixel sampling rates - x direction [px/unit]
    Fs_y = N/Yres;           % Pixel sampling rates - y direction [px/unit]
    
    dx = 1/Fs_x;             % Pixel size - x direction [unit/px]
    dy = 1/Fs_y;             % Pixel size - y direction [unit/px]
    
    dFx = 1/Xres;           % frequency increment - x direction
    dFy = 1/Yres;           % frequency increment - y direction
    
    % Define frequency domain array
            
    Fx = (-M/2:1:M/2-1).*dFx;
    Fy = (-N/2:1:N/2-1).*dFy;

    [GFx, GFy] = meshgrid(Fx,Fy);
    
%% Compute 2D FFT with windowing

    % With windowing
        %{
        w_type = @flattopwin;
        W = window(w_type,N)*window(w_type,M)';
        imgf = fftshift(fft2(W.*Z));    % Fourier transform  
        %}
    % Without windowing
    
        imgf = fftshift(fft2(Z));    % Fourier transform 
        
    
    % Modulus of Fourier coefficients
    
        imgfp = (1/M/N)*abs(imgf);   
    
    % to plot the FFT use:
    % surf(KY,KX,imgfp); view(2); axis equal;
    
                                                                       
%% Compute radial average of 2D FFT


    % For info (we do not consider 2*pi factors here --> to be checked)
    %{
    
    Direct space    |	Frequency space
    -------------------------------
    Pixel size dx   |   Image size 2?/dx
    Pixel size dy   |   Image size 2?/dy
    Image size Lx   |   Pixel size 2?/Lx
    Image size Ly   |   Pixel size 2?/Ly
    Resolution Nx   |   Resolution Nx
    Resolution Mx   |   Resolution Mx
    %}

    if N >= M                     % More rows than columns
        rmax = 0.5/dx;
    elseif N < M                  % More columns than rows
        rmax = 0.5/dy;
    end
    
    
    % Define sampling resolution dr :
    
    dr = hypot(dFx,dFy);
    %dr = 1e4; % fix dr for debug
    r_vec = 0:dr:rmax;
    
    radial_average(1) = imgfp(N/2+1,M/2+1);
    
    for i = 2:numel(r_vec)
        
        radius = r_vec(i);
        num_pts = 2*pi*radius / hypot(dFx,dFy);
        theta = linspace(0,2*pi,num_pts);
    
        xx = radius*cos(theta);
        yy = radius*sin(theta);

        sampled_radial_slice = interp2(GFx,GFy,imgfp,xx,yy);
        radial_average(i) = 2*mean(sampled_radial_slice);  % factor of 2 is because we use single sided fft
        
    end
    
    % For code for debugg
    %{
    
    surf(KY,KX,imgfp); view(2); axis equal; hold on;
    zz = zeros(size(xx));
    plot3(xx,yy,zz,'ro')
    % + change the facet color map to blended (otherwise it seems there is a shift)
    %}

%% Setup output and plot

% Generate plot

figure
plot(r_vec,radial_average,'-ko','LineWidth',1.0)
xlabel('Frequency (au)');
ylabel('Amplitude');
title('Radially averaged Fourier Spectrum');
    set(gca,'XScale','log');
    set(gca,'YScale','log');

% to have a look at the power law

a0 = 120e-6; 
p = 1.7; % power
amps = a0./r_vec.^p;
hold on; plot(r_vec,amps);




