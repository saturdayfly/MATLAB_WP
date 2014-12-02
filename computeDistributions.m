function [z, dz, d2z, p_z, p_dz, p_d2z, cp_dz, cp_d2z, z_nbins, dz_nbins, d2z_nbins] = ...
    computeDistributions(Z, slope_map, curv_map,varargin)
%
% Compute statistical properties of a 3D surface profile
%
% Args: 
%   - Z         : N by M matrice representing surface elevation
%   - slope_map : N by M matrice representing surface slope
%   - curv_map  : N by M matrice representing surface mean curvature
%   - varargin  : optional, a list of property/value pairs
%
% Returns:
%   - z         : elevation values at which probability density p_z has been sampled (based on the range of Z and number of bins z_nbins)
%   - dz        : slope values at which probability density p_dz has been sampled (based on the range of slope_map and number of bins dz_nbins)
%   - d2z       : curvature values at which probability density p_d2z has been sampled (based on the range of curv_map and number of bins d2z_nbins)
%   - p_z       : elevation distribution
%   - p_dz      : slope distribution
%   - p_d2z     : curvature distribution
%   - z_nbins   : number of bins used for discretization of elevation z
%   - dz_nbins  : number of bins used for discretization of slope dz
%   - d2z_nbins : number of bins used for discretization of curvature d2z
%

% Parse arguments
    % Set default values
    options = struct('weights','FALSE');

    optionNames = fieldnames(options);
    % Check count
    if mod(length(varargin),2)~=0
        error('getDistributions needs propertyName/propertyValue pairs')
    end
    % Overwrite options with input values
    for pair = reshape(varargin,2,[])   % pair is {propName;propValue}
        inpName = lower(pair{1});        % make case insensitive
        if any(strcmp(inpName,optionNames))
            % If you want you can test for the right class here
            % Also, if you find out that there is an option you keep getting wrong,
            % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
            options.(inpName) = pair{2};
        else
            error('in getDistributions, %s is not a recognized parameter name',inpName)
        end
    end
    
% Get number of bins to use (user input, if 0 we use Scott rule)
    
    prompt = {  'Number of bins for elevation z', ...
                'Number of bins for slope dz',...
                'Number of bins for curvature d2z'};
    dlg_title = 'Select bining for elevation, slope and curvature histogram computation'; num_lines = 1;
    def = {'0','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    z_nbins     = str2double(answer{1});
    dz_nbins    = str2double(answer{2});
    d2z_nbins   = str2double(answer{3});
    
    % if nbins == 0 we use the Scott rule of thumb
    if z_nbins==0
        z_nbins = floor(range(Z(:))/ (3.5*nanstd(Z(:))*length(Z(:))^(-1/3)));
    end
    if dz_nbins==0
        dz_nbins = floor(range(slope_map(:))/ (3.5*nanstd(slope_map(:))*length(slope_map(:))^(-1/3)));
    end
    if d2z_nbins==0
        d2z_nbins = floor(range(curv_map(:))/ (3.5*nanstd(curv_map(:))*length(curv_map(:))^(-1/3)));
    end

    % Define bin edges
    z_binEdges = linspace(min(Z(:)),max(Z(:)),z_nbins+1);
    dz_binEdges = linspace(min(slope_map(:)),max(slope_map(:)),dz_nbins+1);
    d2z_binEdges = linspace(min(curv_map(:)),max(curv_map(:)),d2z_nbins+1);
    
    % Define bin centers
    z = ( z_binEdges(1:end-1) + z_binEdges(2:end) ) ./ 2;
    dz = ( dz_binEdges(1:end-1) + dz_binEdges(2:end) ) ./ 2;
    d2z = ( d2z_binEdges(1:end-1) + d2z_binEdges(2:end) ) ./ 2;

% Compute univariate distributions

    % height distribution z
    p_z = histc(Z(:),[z_binEdges(1:end-1) Inf])'; p_z(end) = [];
    p_z = p_z/trapz(z,p_z); % normalize

    % slope distribution dz
    if strcmp(options.weights,'TRUE') % slopes are weighted with 1/cos(slope angle)
        
        [~, bin] = histc(slope_map(:),[dz_binEdges(1:end-1) Inf]); % get initial histogram
        weights = 1./cos(atan(slope_map(:)));
        p_dz = accumarray(bin(bin~=0), weights(bin~=0));
        p_dz = (p_dz./trapz(dz,p_dz))';
        
    else % no weight
        
        p_dz = histc(slope_map(:),[dz_binEdges(1:end-1) Inf])'; p_dz(end) = [];
        p_dz = p_dz./trapz(dz,p_dz);
        
    end

    % curvature distribution d2z
    p_d2z = histc(curv_map(:),[d2z_binEdges(1:end-1) Inf])'; p_d2z(end) = [];
    p_d2z = p_d2z./trapz(d2z,p_d2z);
    
    
% Compute bivariate distributions (conditional probabilities)
   
    % slope conditional probability cp_dz
    if strcmp(options.weights,'TRUE') % slopes are weighted with 1/cos(slope angle)
        
        cp_dz = NaN(length(dz),length(z));
        for i=1:length(z_binEdges)-1 % go through each bin
            id = z_binEdges(i) <= Z(:) & Z(:) < z_binEdges(i+1);
            temp_slope  = slope_map(id);
            temp_weight = weights(id);
            [~, bin] = histc(temp_slope(:),[dz_binEdges(1:end-1) Inf]); % get initial histogram
            temp_cp_dz = accumarray(bin(bin~=0), temp_weight(bin~=0),[length(dz) 1]);
            temp_cp_dz = temp_cp_dz./trapz(dz,temp_cp_dz);
            cp_dz(:,i) = temp_cp_dz;
        end
        
    else % no weight
        
        prob = hist3([Z(:) slope_map(:)],'edges',{z_binEdges dz_binEdges});
        prob(end, :) = []; prob(:, end) = [];
        prob = prob'; % for consistency with other methods
        norm = repmat(trapz(dz,prob,1),size(prob,1),1);
        cp_dz =  prob./norm;
        
    end
   
    % curvature conditional probability cp_d2z
    prob = hist3([Z(:) curv_map(:)],'edges',{z_binEdges d2z_binEdges});
    prob(end, :) = []; prob(:, end) = [];
    prob = prob';
    norm = repmat(trapz(d2z,prob,1),size(prob,1),1);
    cp_d2z =  prob./norm;
    
    
    