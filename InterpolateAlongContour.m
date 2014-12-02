function [dz_z, dz2_z, d2z_z] = InterpolateAlongContour(X,Y,Zdata,slope_data,curvature_data,z)
%
% Computes average slope, square slope and curvature as a function of surface
% elevation. Interpolated and average the different quantities along contour lines.
%
% Args:
%   - X, Y              : M by N matrices representing X and Y coordinates in the form of meshgrid
%   - Zdata             : Cell array of M by N matrice containing elevation values
%   - slope_data        : Cell array of M by N matrice containing elevation values
%   - curvature_data    : Cell array of M by N matrice containing elevation values
%   - z                 : vector of elevation values at which quantities are to be interpolated and averaged
%
% Returns:
%   - dz_z      : Average slope function of elevation z
%   - dz2_z     : Average square slope function of elevation z
%   - d2z_z     : Average mean curvature function of elevation z
%

num_dataset = length(Zdata);

% Allocate
Cmat               = cell(1,num_dataset);  % cell array of contour matrices
Cdat               = cell(1,num_dataset);  % cell array of contour data ( = Cmat arranged in structures)
S                  = cell(1,num_dataset);  % cell array of slope and curvature data
dz_z_temp    = NaN(length(z),num_dataset); % average slope dz function of z (columns corresponds to different dataset)
dz2_z_temp   = NaN(length(z),num_dataset); % average square slope dz2 function of z (idem)
d2z_z_temp   = NaN(length(z),num_dataset); % average curvature d2z function of z (idem)

% Extract coordinates of contour lines
for i=1:num_dataset
    fprintf('Computing contour for dataset number %d ...\n',i);
    Cmat{i} = contourc(X(1,:),Y(:,1),Zdata{i},z);     % compute contour
    Cdat{i} = reshapeContourData(Cmat{i},z);          % reshape data
end

% Intepolate the slope & curvature along the contour lines
for i=1:num_dataset
    for j = 1:length(z)
        
        S{i}(j).dz 	= interp2(X,Y,slope_data{i},Cdat{i}(j).xdata,Cdat{i}(j).ydata,'linear');
        S{i}(j).dz2 = S{i}(j).dz.^2;
        S{i}(j).d2z = interp2(X,Y,curvature_data{i},Cdat{i}(j).xdata,Cdat{i}(j).ydata,'linear');
        
        S{i}(j).dz_mean     = nanmean(S{i}(j).dz);
        S{i}(j).dz_std      = nanstd(S{i}(j).dz);
        S{i}(j).dz2_mean    = nanmean(S{i}(j).dz2);
        S{i}(j).dz2_std     = nanstd(S{i}(j).dz2);
        S{i}(j).d2z_mean    = nanmean(S{i}(j).d2z);
        S{i}(j).d2z_std     = nanstd(S{i}(j).d2z);
        
        dz_z_temp(j,i)    = S{i}(j).dz_mean;
        dz2_z_temp(j,i)   = S{i}(j).dz2_mean;
        d2z_z_temp(j,i)   = S{i}(j).d2z_mean;
        
    end
end

% Compute mean slope, slope² and curvature by combining the datasets
dz_z  = nanmean(dz_z_temp,2)';
dz2_z = nanmean(dz2_z_temp,2)';
d2z_z = nanmean(d2z_z_temp,2)';

