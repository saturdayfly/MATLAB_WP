function [x0, y0, dy] = findJump(x,y,xres,dx,lower_bound,upper_bound)

% find discontinuity along a curve y = f(x)
% used to identify jumps in adsorbed liquid in the desorption isotherms
%
% Args:
%   - x : x coordinate of the curve (typically the curvature)
%   - y : y coordinate of the curve (typically average film position)
%   - xres : resolution of the x data. Used to compute the amplitude of a
%            jump
%   - dx : used to linearly resample the x data (x is potentially a log spaced variable)
%   - lower_bound,upper_bound: minium and maximum y values where a jump is
%                              to be looked for
% Returns:
%   - x0 : x position of the discontinuity
%   - y0 : y position of the discontinuity
%   - dy : amplitude of the discontinuity =  abs(y(x0+xres/2) - y(x0-xres/2))

% possible optimization : work directly with indices as x coordinate.

% spline interpolate the curve and find the strongest inflection point
% with negative derivative

    % spline interpolate
    xx = min(x):dx:max(x);
    cs = spline(x,y);

    % derivatives
    fprime  = fnder(cs,1);
    fsecond = fnder(cs,2);

    % find all inflection points within bounding box
    inflec_pts = fnzeros(fsecond);
    
    inflec_pts(:,inflec_pts(1,:)~=inflec_pts(2,:)) = []; % remove zero splines
    inflec_pts(2,:) = [];                                % remove second raw
    inflec_pts(ppval(cs,inflec_pts) > upper_bound) = [];   % remove inflection points corresponding to data lying above the upper bound
    inflec_pts(ppval(cs,inflec_pts) < lower_bound) = [];   % remove inflection points corresponding to data lying below the lower bound

    % look for the inflection point with strongest negative derivative
    if ~isempty(inflec_pts) % found a discontinuity
            
        derivatives = ppval(fprime,inflec_pts); % compute corresponding derivative
        [~, id2] = min(derivatives); % pick indice corresponding to largest negative slope
        
        x0 = inflec_pts(id2); % x value at which the jump occurs
        y0 = ppval(cs,x0);
               
        % y value on the left side of the jump
        [~, idup] = min(abs(x-(x0-xres/2)));
        y1 = y(idup);
            
        % y value on the right side of the jump
        [~, iddown] = min(abs(x-(x0+xres/2)));
        y2 = y(iddown);
        
        % amplitude of the jump
        dy = y1-y2;
        
    else % did not find any discontinuity, returns NaN
   
        x0 = NaN;
        y0 = NaN;
        dy = NaN;
        
    end
