%	Script to compute Wetting phase diagram for a randomly rough surface
%	Written by R. DUFOUR - MPIDS/DMI - renaud.dufour@ds.mpg.de
%	Version : Decmeber 2014
%
%   This script is based on the following papers :
%
%   [1] S. Herminghauss, PRL, 109, 2012
%       "Universal Phase Diagram for Wetting on Mesoscale Roughness"
%
%   [2] S. Herminghauss, Eur. Phys. J. E., 35(43), 2012
%       "Wetting, Spreading and adsorption on randomly rough surfaces"
%
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   - PART 1    : Load or generate numerically Z profiles
%
%   - PART 2.1  : compute surface slope and curvature
%   - PART 2.2  : compute surface statistical properties
%   - PART 2.3  : compute sigma functions (i.e. dz,dz2 & d2z vs z)
%
%   - PART 3.1   : Compute desorption isotherms
%   - PART 3.2   : Compute phase diagram from desorption isotherms
%
%   - PART 4     : Data export

%% Setup

format long;
addpath('surface generators/')

%% Part 1 : load or generate surface profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load from file
%[X, Y, Zdata, surface_attributes] = importFromGwyddion();

% generate - wrapper
%{
[X, Y, Zdata{1}, attributes] = makeSurface('rough',...
                                        'num_modes',5,...
                                        'num_points',100,...
                                        'power',2);
%}

% generate - custom
%for this first load the phase data
%{
[X, Y, Zdata{1}, attributes] = generateRough20modes(20,400,phase);
num_dataset = numel(Zdata);
%}

% To plot use : surf(X,Y,Zdata);

%% Part 2.1 : compute slope and curvature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Compute Slope and curvature...');

% parameters
ROI_size     = 3;            % size of moving window for computeSlope()
slope_method = 'mldivide';   % method for computeSlope()

% allocate
slope_data     = cell(1,num_dataset);
curvature_data = cell(1,num_dataset);

% compute slope & curvature
for k=1:num_dataset
    slope_data{k}           = computeSlope(X,Y,Zdata{k},ROI_size,slope_method);
    [~ , curvature_data{k}] = computeCurvature(X,Y,Zdata{k});
end

% concatenate the data
Z_map     = cell2mat(Zdata);
slope_map = cell2mat(slope_data);
curv_map  = cell2mat(curvature_data);

fprintf('Done.\n');
    
%% Part 2.2 : compute surface statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Compute Statistical properties...');

% Compute some global statistical properties
rms_amplitude   = sqrt(var(Z_map(:)));      % Root mean square amplitude
rms_slope       = nanstd(slope_map(:));     % Root mean square slope
rms_curvature   = sqrt(var(curv_map(:)));   % Root mean square curvature

skew        = skewness(Z_map(:));           % skewness
kurt        = kurtosis(Z_map(:))-3;         % normalized kurtosis

% Add correlation length




% Select a method to compute statistical quantities
%   - histogram : histogram based approach
%   - weighted :  same than histogram, but slopes are weighted with 1/cos
stat_method = 'histogram';
% Activate weighting of slope (implemented for histogram method only)
use_weights = 'FALSE';

switch stat_method
    case 'histogram'
        [z, dz, d2z, p_z, p_dz, p_d2z, cp_dz, cp_d2z, z_nbins, dz_nbins, d2z_nbins] = ...
            computeDistributions(Z_map, slope_map, curv_map,'weights',use_weights);
    case 'kernel'
        [z, dz, d2z, p_z, p_dz, p_d2z, cp_dz, cp_d2z, global_nbins] = ...
            kdens_stat(Z_map, slope_map, curv_map);
end

figure;
subplot(1,3,1); plot(z,p_z,'-bo'); % Height distribution
title('Height distribution (normalized)'); xlabel('h');
subplot(1,3,2); plot(dz,p_dz,'-bo'); % Slope distribution
legend(sprintf('Plane fit (%d*%d pixels)',ROI_size,ROI_size))
title('Slope distribution (normalized)'); xlabel('Slope');
subplot(1,3,3); plot(d2z,p_d2z,'-ro'); % Curvature distribution
legend(sprintf('Plane fit (%d*%d pixels)',ROI_size,ROI_size))
title('Curvature distribution (normalized)'); xlabel('Curvature');

fprintf('Done.\n');
    
%% Part 2.3 : compute average slope & curvature as a function of height

fprintf('Compute Sigma functions...');

% Select a method to compute the sigma functions
%    - bining: rely on the slope conditional probability
%    - interpolate: interpolate and average slopes along contour lines
sigma_method = 'bining';

switch sigma_method
    case 'bining'
        dz_z    = NaN(size(z));     % average slope dz function of z
        dz2_z   = NaN(size(z));     % average square slope dz2 function of z
        d2z_z   = NaN(size(z));     % average curvature d2z function of z
        for h=1:1:length(z)
            dz_z(h)     = trapz(dz,cp_dz(:,h)'.*dz);
            dz2_z(h)    = trapz(dz,cp_dz(:,h)'.*dz.^2);
            d2z_z(h)    = trapz(d2z,cp_d2z(:,h)'.*d2z);
        end
    case 'interpolate'
        [dz_z, dz2_z, d2z_z] = InterpolateAlongContour(X,Y,Zdata,slope_data,curvature_data, z);
end

figure('name','Average slope, square slope and curvature function of elevation','units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1);
plot(z, dz_z,'-ob',z,dz2_z,'-or')
legend('Average slope vs elevation (sigma1)','Average square slope vs elevation (sigma2)');
xlabel('z'); ylabel('sigma1, sigma2');

subplot(1,3,2);
plot(z,dz2_z./dz_z,'-ok')
legend('sigma2/sigma1'); xlabel('z'); ylabel('sigma2/sigma1');

subplot(1,3,3);
plot(z(d2z_z>0),d2z_z(d2z_z>0),'-ob'); hold on;
plot(z(d2z_z<0),abs(d2z_z(d2z_z<0)),'-or'); hold off;
legend('Positive curvature','Negative curvature'); xlabel('z'); ylabel('Curvature');
set(gca,'YScale','log');

fprintf('Done\n');

%% part 3.1 : Compute desorption isotherms

fprintf('Compute Sigma functions...');

% User input for data processing
prompt = {  'Threshold for phase diagram computation [%]', ...  % Determine the theshold below which the surface is considered as dry ...
            'Theta increment in degree(default 1)',...          %   ... it can be used to compute percolation diagram setting threshold = 0 for example.
            'Curvature : number of values (default 500)',...
            'Minimal curvature (in power of 10)',...
            'Display graphical resolution (on/off)',...
            'Maximum curvature (default : max curv of surface)',...
            'method for smoothing',...
            'span for smoothing'};

dlg_title = 'Analasis parameters'; num_lines = 1;
def = {'0','1','500','-5','off',num2str(max(abs(curv_map(:)))),'lowess',num2str(0.1)};

answer = inputdlg(prompt,dlg_title,num_lines,def);
threshold        = str2double(answer{1});
theta_inc        = str2double(answer{2});
H_nvalues        = str2double(answer{3});
p                = str2double(answer{4});
display          = answer{5};
H_max            = str2double(answer{6});
smooth_method    = answer{7};
smooth_span      = str2double(answer{8});

% Define functions sigma1, sigma2, L and W
sigma1 = dz_z;               % momentum of order 2 of dz
sigma2 = dz2_z;              % momentum of order 2 of dz2
L = sigma1.*p_z;             % length of contour line a height z
W = NaN(size(z)); W(1) = 0;  % wetted sample area
for h=2:1:length(z)
    W(h) = trapz(z(1:h),p_z(1:h));
end

% Compute Wenzel angle
thetaw = sqrt(trapz(z,sigma2.*p_z));

% Solve
data1 = computeIsotherms( p_z,sigma1,sigma2,L,W,z,...
              H_max,theta_inc,H_nvalues,p,...
              threshold,smooth_method,smooth_span,display);

% Plot Results
np = 2; % plot every 'np' value of theta

% Average film position vs TOTAL curvature 2*H
figure('name','Average film position (METHOD 1)'); hold on;
plot(2*data1.H_vec,data1.AVG_FILM_POSITION(1:np:end,:),'marker','o');
count = 1;
for t = 1:np:length(data1.theta_vec) % define legend
    graphlegend{count} = sprintf('\\theta = %1.1f',atand(data1.theta_vec(t)));
    count = count + 1;
end
legend(graphlegend); xlabel('Curvature H'); ylabel('Average film position');
plot([0 2*max(data1.H_vec)],[threshold threshold],'-'); % percolation threshold
hold off;
set(gca,'XScale','log');

% Adsorbed volume per unit area vs TOTAL curvature 2*H
figure('name','Adsorbed volume (METHOD 1)'); hold on;
plot(2*data1.H_vec,data1.VOLUME(1:np:end,:),'marker','o');
legend(graphlegend); xlabel('Curvature H'); ylabel('Adsorbed volume');
hold off;
set(gca,'XScale','log');
%set(gca,'YScale','log');

% Percolation diagram (angle vs TOTAL curvature 2*H)
figure('name','Percolation diagram (METHOD 1)'); hold on;
plot(2*data1.Hp,atand(data1.theta_vec), 'ob');
plot(0,atand(data1.thetaf),'or',0,atand(thetaw),'or');
title('Percolation diagram'); xlabel('Curvature H'); ylabel(sprintf('\\theta [Deg]'));
text(0.001,atand(thetaw),'\theta_W','FontSize',12)
hold off;

% To plot a graphical resolution for a given angle and curvature use :
% plotSolution(sigma1,sigma2,L,W,z,10,0.02) 

fprintf('Done\n');

%% part 3.2 : Compute phase diagram  

% parameters
h_error     = 1e-5;
dh          = 1e-6;

phase_diagram = computePhaseDiagram(data1,h_error,dh);

% plot phase diagram
figure('name','Phase diagram');
plot(2*phase_diagram.curvature,atand(phase_diagram.theta),'o');
    

%% part 4 : Export data


    % To to : put sample name and date in the filename.
    prompt = {'Sample name'};
    dlg_title = 'Input sample name for data export';
    def = {'noname'};
    answer = inputdlg(prompt,dlg_title,1,def);
    prefix = sprintf('%s_%s_',datestr(now,29),answer{1});

    % delete dat files
    delete( '*.dat');
    
    dest = 'results';
    if ~isdir(dest)
        mkdir(dest)
    end
    
    % Export
    Script_exportData

    
% end