function [z, dz, d2z, p_z, p_dz, p_d2z, cp_dz, cp_d2z, global_nbins] = ...
    kdens_stat(Z, slope_map, curv_map)

% Compute statistical properties of a 3D surface profile based on kernel
% density estimators
%
% Args: 
%   - Z         : N by M matrice representing surface elevation
%   - slope_map : N by M matrice representing surface slope
%   - curv_map  : N by M matrice representing surface mean curvature
%
% Returns:   
%   - z             : elevation values at which probability density p_z has been sampled
%   - dz            : slope values at which probability density p_dz has been sampled
%   - d2z           : curvature values at which probability density p_d2z has been sampled
%   - p_z           : pdf of elevation
%   - p_dz          : pdf of slope
%   - p_d2z         : pdf of curvature
%   - cp_dz         : conditional pdf of slope 
%   - cp_d2z        : conditional pdf of curvature
%   - global_nbins  : number of bins used for discretization of z, dz and d2z
%

    prompt = {  'Number of bins (power of 2)'};
    dlg_title = 'Select bining for kernel density estimation'; num_lines = 1;
    def = {'256'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    global_nbins     = str2double(answer{1});
    
    % Compute joint probability of height & slope using 2d kernel density
    data_z_dz = [Z(:) slope_map(:)];
    MIN_XY = [min(Z(:)) 0];
    MAX_XY = [max(Z(:)), max(slope_map(:))];
    [~,density_z_dz,mesh_z,mesh_dz] = kde2d(   data_z_dz, ...
                                                            global_nbins, ...
                                                            MIN_XY, ...
                                                            MAX_XY);
                                                        
    % Extract height and slope probabilities
    z    = mesh_z(1,:);
    dz   = mesh_dz(:,1)';
    p_z  = trapz(dz,density_z_dz,1);
    p_dz = trapz(z,density_z_dz,2)';
    
    % Extract conditional slope probabilities + plot
    cp_dz =  density_z_dz./repmat(p_z,size(density_z_dz,1),1);
    
    % Compute joint probability of height & curvature using 2d kernel density
    data_z_d2z = [Z(:) curv_map(:)];
    MIN_XY = [min(Z(:)) min(curv_map(:))];
    MAX_XY = [max(Z(:)), max(curv_map(:))];
    [~,density_z_d2z,mesh_zbis,mesh_d2z] = kde2d( data_z_d2z, ...
                                                                global_nbins, ...
                                                                MIN_XY, ...
                                                                MAX_XY);
    
    % Extract height and slope probabilities
    zbis    = mesh_zbis(1,:);
    d2z     = mesh_d2z(:,1)';
    p_zbis  = trapz(d2z,density_z_d2z,1);
    p_d2z    = trapz(zbis,density_z_d2z,2)';
    
    % Extract conditional curvature probabilities + plot
    cp_d2z =  density_z_d2z./repmat(p_zbis,size(density_z_d2z,1),1);
    
    