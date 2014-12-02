function S = reshapeContourData(c,z)
% Reshape the contour data returned by matlab contour functions.
% CONTOUR, CONTOURF, CONTOUR3, and CONTOURC all produce a contour matrix
% C that is traditionally used by CLABEL for creating contour labels.

% The contour matrix C is a two row matrix of contour lines.
% Each contiguous drawing segment contains the value of the contour,
% the number of (x, y) drawing pairs, and the pairs themselves.
% The segments are appended end-to-end as :
 
%        C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
%             pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]
             
%
% S = reshapeContourData(c,z) extracts the (x,y) data pairs describing each
% contour level and other data from the contour matrix C.
%
% Args: 
%   - c: the contour matrix returned by CONTOUR CONTOURF CONTOUR3 or CONTOURC
%   - z: vector of elevation values at which contour lines were computed
%     (providing z in addition to c enables handling of empty contours)
%
% Returns:
%   - S: a structure array of length length(z) with the following fiels :
%       S(k).level: the contour level height. Equals z(k)
%       S(k).numel: number of points describing the level k.
%       S(k).xdata: x-axis data for the level k.
%       S(k).ydata: y-axis data for the level k.
%
% For example: PLOT(S(k).xdata,S(k).ydata)) plots all contours at level k.
%
% Notes :
%
% This code was inspired by the CONTOURDATA function of
% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
%
% The main difference is that all the (x,y) values belonging to the same
% elevation are concatenated in a single S(i).
% Also the list of contour levels z is used as input so that we can ensure
% that size(S) = size(z) and have a correspondance z(i) == S(i)
% (otherwise, there would be no structure S(i) returned for 'empty' levels)
%
% Improvements (todo) :
%
% 	-as a preprocessing step, interpolate the data along each contour
%    to get evenly spaced points. Use for example interparc function :
%    tt = interparc(length(x),x,y); % default is spline fitting

if nargin<1 || ~isfloat(c) || size(c,1)~=2 || size(c,2)<4
   error('CONTOURDATA:rhs',...
         'Input Must be the 2-by-N Contour Matrix C.')
end

% initialize structure s
S = struct( 'level',num2cell(z), ...
            'numel',num2cell(zeros(size(z))),...
            'xdata',[],...
            'ydata',[]);

col = 1;        % index of column containing contour level and number of points

for i = 1:length(z) % go through each level
    
    while col<size(c,2) && c(1,col) == z(i);    % while we find contour data corresponding to level i
                                                % and we have not reached the end of C
        
        % TO DO : preprocessing step to get evenly spaced points along
        % individual contours
                          
        S(i).numel  = S(i).numel + c(2,col);
        idx         = col+1:col+c(2,col);
        S(i).xdata  = [S(i).xdata; c(1,idx).'];   % concatenate
        S(i).ydata  = [S(i).ydata; c(2,idx).'];   % concatenate
        
        col=col+c(2,col)+1; % go to next contour
    end
    
    % we go out of the while loop, means we jump to the next level
    
end



