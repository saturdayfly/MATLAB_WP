function [X, Y, Z, attributes] = makeSurface(type,varargin)
%
% Generate a surface profile numerically (wrapper)
%
% Args:
%   - type: the type of profile to be generated. types supported are :
%           "rough", "egg cartoon", "egg cartoon hexa", "hole", "hole 2D"
%           "wave", "multimode 2D"
%   - a list of property/value pairs specific to the type of function
%   called. See table below:
%


% Parse arguments
    % Set default values
    options = struct(   'num_points',NaN,...
                        'num_modes',NaN,...
                        'power',NaN,...
                        'max_slope',NaN,...
                        'rep',NaN);

    optionNames = fieldnames(options);
    % Check count
    if mod(length(varargin),2)~=0 
       error('makeSurface needs propertyName/propertyValue pairs')
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
          error('in makeSurface, %s is not a recognized parameter name',inpName)
       end
    end

% Call subfunction according to type
switch type
    case 'rough'
        if (isnan(options.num_points) || isnan(options.num_modes) || isnan(options.power))
            error('in makeSurface, invalid options for type %s.\n',type)
        else
        [X, Y, Z] = generateRough(  options.num_modes,...
                                    options.num_points,...
                                    options.power);
        end
    case 'egg cartoon'
        if (isnan(options.num_points) || isnan(options.max_slope))
            error('in makeSurface, invalid options for type %s.\n',type)
        else
        [X, Y, Z] = generateEggCartoon(options.max_slope,...
                                       options.num_points);
        end
    case 'egg cartoon hexa'
        if (isnan(options.num_points) || isnan(options.max_slope))
            error('in makeSurface, invalid options for type %s.\n',type)
        else
        [X, Y, Z] = generateEggCartoonHexa( options.max_slope,...
                                            options.num_points);
        end
    case 'hole'
        if (isnan(options.num_points) || isnan(options.max_slope))
            error('in makeSurface, invalid options for type %s.\n',type)
        else
        [X, Y, Z] = generateHole( options.max_slope,...
                                  options.num_points);
        end
    case 'hole 2D'
        if (isnan(options.num_points) || isnan(options.max_slope) || isnan(options.rep))
            error('in makeSurface, invalid options for type %s.\n',type)
        else
        [X, Y, Z] = generateHole2D( options.max_slope,...
                                    options.num_points,...
                                    options.rep);
        end
    case 'wave'
        [X, Y, Z] = generateWave();
    case 'multimode 2D'
        if (isnan(options.num_points) || isnan(options.max_slope) || isnan(options.rep))
            error('in makeSurface, invalid options for type %s.\n',type)
        else
        [X, Y, Z]    = generateRough2D(options.num_modes,...
                                       options.num_points,...
                                       options.rep);
        end
    otherwise
        error('in makeSurface, %s is not a recognized type',type)
end

    attributes.width        = 1;
    attributes.height       = 1;
    attributes.x_resolution = 1/options.num_points;
    attributes.y_resolution = 1/options.num_points;
    attributes.unit         = '1';