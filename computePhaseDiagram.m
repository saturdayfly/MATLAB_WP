function [phase_diagram] = computePhaseDiagram(data,h_error,dh)
% Compute phase diagram from desorption isotherms
%
% Args:
%   - data: dataset obtained from WPD1 function
%   - h_error: assumed error on the curvature (passed to findJump)
%   - dh: curvature discretization step (passed to findJump)
%
% returns:
%   - phase_diagram: structure containing the following fields:
%       .curvature: curvature at which jump occurs
%       .position: film position at which jump occurs
%       .amplitude: amplitude of the jump
%       .theta: contact angle at which the jump occurs

% Note : phase diagram is computed from isotherms representing
% pressure vs average film position. Thus phase_diagram.position represents
% an average position of the liquid film, and .amplitude unit is a distance.
% As an alternative, one may want to implement the alternative computation
% considering desorption isotherms in the form of pressure vs adsorbed volume.

% parameters
angle     = data.theta_vec; % vector of contact angle
curvature = data.H_vec;     % vector of curvature values (mean curvature)

% We first limit the observation windows to make computation easier

    % Select window (user input)
    fig = figure;
    plot(curvature,data.AVG_FILM_POSITION(:,:),'marker','o');
    set(gca,'XScale','log');

    dcm_obj = datacursormode(fig);
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
    disp('Select 2 points (hold Alt key), then press Return.')
    pause;
    c_info = getCursorInfo(dcm_obj); % contains cursor infos
    close(fig);

    % retrieve cursor information
    if c_info(1).Position(1) < c_info(1).Position(1)
        cursor1 = c_info(1);
        cursor2 = c_info(2);
    else
        cursor1 = c_info(2);
        cursor2 = c_info(1);
    end

    % Define bounding box
    idx1 = cursor1.DataIndex;
    idx2 = cursor2.DataIndex;
    lower_bound = min([cursor1.Position(2) cursor2.Position(2)]);
    upper_bound = max([cursor1.Position(2) cursor2.Position(2)]);

% Go through each desorption isotherm and look for jump in film thickness

    % Allocate
    pd.curvature     = nan(length(angle),1); % curvature at which jump occur
    pd.position      = nan(length(angle),1); % film position at which the jump occurs
    pd.amplitude     = nan(length(angle),1); % amplitude of the jump
    
    % we display the solution while solving
    jumpCalculation = figure('name','jump calculation');
    
    xdata = curvature(idx1:idx2);
    for j = 1:length(angle)
       
        ydata = data.AVG_FILM_POSITION(j,idx1:idx2);
        [pd.curvature(j), pd.position(j), pd.amplitude(j)] = ...
            findJump(xdata,ydata,h_error,dh,lower_bound,upper_bound);
        
        % live plot
        hold on;
        plot(xdata,ydata,'-b');
        plot(pd.curvature(j),pd.position(j),'Marker','p','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
        plot([pd.curvature(j)-h_error/2 pd.curvature(j)+h_error/2],[pd.position(j) pd.position(j)],'-xk');
        plot([pd.curvature(j) pd.curvature(j)],[pd.position(j)-pd.amplitude(j)/2 pd.position(j)+pd.amplitude(j)/2],'-xk');
        set(gca,'XScale','log'); axis tight;
        hold off;
        pause(0.3);
        drawnow;
    end

    close(jumpCalculation);
    phase_diagram = pd;
    phase_diagram.theta = angle;
    
end

