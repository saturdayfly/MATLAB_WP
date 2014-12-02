function [X, Y, Zdata, attributes] = importFromGwyddion()
%
% Import surface profile from file. 
% We assume tab delimited text files with no header.
% Note: from Gwyddion, export in Ascii format without header.
%
% Returns :
%   - X, Y: coordinates, two N by M matrices created with meshgrid
%   - Zdata: a cell array containing N by M matrices with surface elevation
%            length of Zdata correspond to the number of input files   
%   - attributes : a structure containing the following variables :
%       .unit: the unit of X,Y and Z coordinates
%       .width: width of the surface in unit (x dimension)
%       .height: height of the surface in unit (y dimension)
%       .x_resolution: resolution in unit per pixel
%       .y_resolution: resolution in unit per pixel

delimiter = sprintf('\t'); % tab delimited files

% Select files      
[filename, pathname] = uigetfile({'*.*'},'Select files','MultiSelect','on');
if ischar(filename);
    filelist{1} = filename; % coerce char to cell
else
    filelist = filename;
end
ndataset = numel(filelist);

% Ask User file informations
prompt = {  'Image width [m]', ... % Image width in pixels
            'Image height [m]',... % Image height in pixels
            'Target unit',...      % Unit of data (after conversion factor is applied)
            'Conversion factor'};  % conversion factor to match unit
dlg_title = 'Files properties'; num_lines = 1;
def = {'600e-6','450e-6','µm','1e6'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

rescale_factor  = str2num(answer{4});
width           = str2num(answer{1})*rescale_factor;
height          = str2num(answer{2})*rescale_factor;
unit            = answer{3};

% Determine number of columns
fid         = fopen(fullfile(pathname,filelist{1}),'rt');
first_line  = fgets(fid);
num_columns = numel(strfind(first_line,delimiter)) + 1;
fclose(fid);

% Read files
formatSpec  = repmat('%f',1,num_columns); % note: in matlab generated import 
                                          %       function, formatSPec ends
                                          %       with [^\n\r]  (?)
Zdata       = cell(1,ndataset);
for k=1:ndataset
    current_file = fullfile(pathname,filelist{k});
    fileID = fopen(current_file,'r');
    dataArray    = textscan( fileID, formatSpec, inf,...
                            'Delimiter', delimiter,...
                            'EmptyValue' ,NaN,...
                            'HeaderLines', 0,...
                            'ReturnOnError', false);
    Zdata{k}    = [dataArray{1:end}];
    Zdata{k}    = Zdata{k}.*rescale_factor; % convert to target unit
    fclose(fileID);
end
num_rows = size(Zdata{1},1);

% generate grid
dx = width/num_columns;
dy = height/num_rows;
[X,Y] = meshgrid(dx*(1:1:num_columns), dy*(1:1:num_rows));
    
attributes.width        = width;
attributes.height       = height;
attributes.x_resolution = dx;
attributes.y_resolution = dy;
attributes.unit         = unit;
    
% To plot imported surface profile use : surface(X,Y,Zdata{k}); view(3);
    
% End       
    
        
% set pixel size size according to magnification (see Weeko NT1100 doc)
% From previous version, not used anymore
%{
    if strcmp(mag,'10')
        dx      = 0.815*1e-6;     % pixelsize in [m] - x direction
        dy      = 0.937*1e-6;     % pixelsize in [m] - y direction
    elseif strcmp(mag,'2.5')
        dx      = 3.329*1e-6;
        dy      = 3.986*1e-6;
    elseif strcmp(mag,'0')        % I use this for surface profiles generated with gwyddion
        dx      = 0.815*1e-6;     % specify below the pixel size used (usually square pixels)
        dy      = 0.815*1e-6;
    else
        fprintf('Error: Unknown/Uncorrect magnification. Select 2.5 or 10.\n'); break;
    end
%}


    