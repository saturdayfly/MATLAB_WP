function [R, rho] = computeRadialACF(X,Y,Z,drho,scaleType)
%
% Compute the radial autocorrelation function of a 3D surface.
% STILL HAVE TO DOUBLE CHECK THE OUTPUT COMPARING WITH GWYDDION
%
% Args: 
%   - X, Y : M by N matrices, x and y coordinates (from meshgrids)
%   - Z : M by N matrix : z coordinate
%   - scaleType: character string, type of normalization
%       biased   : scales the raw auto-correlation by the maximum number
%                  of elements involved its generation
%       unbiased : scales the raw auto-correlation by dividing each element
%                  by the number of products A and B used to generate that
%                  element
%       coeff    : normalizes the sequence so that the largest
%                  auto-correlation element equals 1.0
%       none     : no scaling
%
% Returns :
%   - R    : row vector containing the autocorrelation data
%   - rho  : row vector containing the lag values
%

[nrow,ncol] = size(Z);

% Compute 2D autocorrelation of Z with xcorr2 function
% note: output is not normalized - see help for xcorr2

    hi = popup; watchon; drawnow;   % popup window during computation
    C = xcorr2(Z);
    watchoff; close(hi);            % close popup

% Normalise the autocorrelation

    switch scaleType,
        case 'none',
            % do nothing
        case'biased',
            % Scales the raw cross-correlation by 1/(M*N)
            C = C/(nrow*ncol);
        case 'unbiased',
            % scales the raw auto-correlation by dividing each element by
            % the number of products A and B used to generate that element
            [M,N] = size(Z);
            C = C./([1:M-1 M:-1:1]'*[1:N-1 N:-1:1]);
            % see www.mathworks.com/matlabcentral/newsreader/view_thread/59030
        case'coeff',
        % Normalizes the sequence so that the auto-correlations
        % at zero lag are identically 1.0.
            C = C./max(C(:));
    end

% Average the 4 quadrants (assume isotropic profile):

    c1 = C(size(Z,1):end,size(Z,2):end);
    c2 = C(size(Z,1):end,1:size(Z,2));      c2 = fliplr(c2);
    c3 = C(1:size(Z,1),1:size(Z,2));        c3 = rot90(c3,2); %180° rotation
    c4 = C(1:size(Z,1),size(Z,2):end);      c4 = flipud(c4);
    Zcorr_cart = (c1+c2+c3+c4)/4; % averaging

% Define lag grids and convert to polar coordinates

    dx = X(1);
    dy = Y(1);
    [XLAG, YLAG]            = meshgrid(0:dx:(ncol-1)*dx, 0:dy:(nrow-1)*dy);
    [THETA,RHO,Zcorr_polar] = cart2pol(XLAG,YLAG,Zcorr_cart);

% Interpolate onto a regular polar grids

    rho_max = min([ max(X(:)) max(Y(:)) ]);
    dtheta  = drho/rho_max;
    [THETA_int, RHO_int] = meshgrid( 0:dtheta:max(THETA(:)), ...
                                     0:drho:rho_max);

    Zcorr_polar_int = griddata(THETA(:),RHO(:),Zcorr_polar(:),THETA_int,RHO_int);

% Finally, compute the angular average

    R       = nanmean(Zcorr_polar_int,2);
    rho     = nanmean(RHO_int,2);

    function [hi] = popup
        
        sz = [250 100]; % figure size
        screensize = get(0,'ScreenSize');
        xpos = ceil((screensize(3)-sz(2))/2);
        ypos = ceil((screensize(4)-sz(1))/2);
        hi = figure(...
            'position',[xpos, ypos, sz(1), sz(2)],...
            'units','pixels',...
            'renderer','OpenGL',...
            'MenuBar','none',...
            'PaperPositionMode','auto',...
            'color','black',...
            'Name','Data Input Dialog',...
            'NumberTitle','off',...
            'Tag','gui',...
            'Resize','off');
        
        h = uicontrol( ...
            'Style','text', ...
            'units','pixels',...
            'position',[0 30 250 30],...
            'FontSize',13,...
            'String','Computing Auto-correlation ...',...
            'BackgroundColor','black',...
            'ForegroundColor','white');
        