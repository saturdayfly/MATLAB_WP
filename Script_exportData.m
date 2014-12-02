% A few shortcuts to make things simpler
    unit = attributes.unit;
    % Export parameters

    filename = fullfile(dest,sprintf('%s%s',prefix,'parameters.dat'));
    
    header_longname={   'With of datafile',...
                        'Height of datafile',...
                        'unit',...
                        'lateral pixel size',...
                        'vertical pixel size',...
                        'number of dataset',...
                        'z number of bins',...
                        'slope number of bins',...
                        'curvature number of bins',...
                        'size of fitted plane for slope computation',...
                        'method of slope computation',...
                        'discretization step for theta (linearly spaced)',...
                        'number of curvature values used to compute desorption isotherm (log spaced)',...
                        'maximum mean curvature considered',...
                        'mean curvature error considered for jump computation',...
                        'discretization step for the curvature for jump computation',...
                        'Wenzel angle',...
                        };
        
    header_unit={       'px','px','','',unit,unit,'','','','','px','','deg','',...
                        sprintf('1/%s',unit),sprintf('1/%s',unit),sprintf('1/%s',unit),'rad'}; 
                    
    header_comment={    '','','','','','','','','','','','','','','','','','',''}; 
                    
    parameters = {  num2str(attributes.width),...
                    num2str(attributes.height),...
                    attributes.unit,...
                    num2str(attributes.x_resolution),...
                    num2str(attributes.y_resolution),...
                    num2str(num_dataset),...
                    num2str(length(z)),...
                    num2str(length(dz)),...
                    num2str(length(d2z)),...
                    num2str(ROI_size),...
                    slope_method,...
                    num2str(theta_inc),...
                    num2str(H_nvalues),...
                    num2str(H_max),...
                    num2str(h_error),...
                    num2str(dh)...
                    num2str(thetaw)
                    };
                   
    txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
    txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
    txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';
    txt4=sprintf('%s\t',parameters{:}); txt4(end)='';
    dlmwrite(filename,txt1,'delimiter','');
    dlmwrite(filename,txt2,'delimiter','','-append');
    dlmwrite(filename,txt3,'delimiter','','-append');
    dlmwrite(filename,txt4,'delimiter','','-append');
    
    % Export statistical properties

    % statistics 1 : z, p(z), sigma1, sigma 2, average curvature(absolute) + sign and L
        filename = fullfile(dest,sprintf('%s%s',prefix,'statistic1.dat'));

        header_longname ={ 'z',     'z distribution',       'sigma1',   'sigma2',   'Average curvature (absolute value)',   'Average curvature (sign)', 'Length of contour line (per unit area)'  };
        header_unit     ={ unit,    sprintf('1/%s',unit),   '',         '',         sprintf('1/%s',unit),                   '',                         sprintf('1/%s',unit)};
        header_comment  ={ '',      '',                     '',         '',         '',                                     '',                         ''};

        txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
        txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
        txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';

        dlmwrite(filename,txt1,'delimiter','');
        dlmwrite(filename,txt2,'delimiter','','-append');
        dlmwrite(filename,txt3,'delimiter','','-append');

        dlmwrite(filename,[z' p_z' sigma1' sigma2' d2z_z' sign(d2z_z'), L'],'delimiter','\t','-append');

    % statistics 2 : slope distribution
        filename = fullfile(dest,sprintf('%s%s',prefix,'statistic2.dat'));

        header_longname ={'slope',  'slope distribution'    };
        header_unit     ={'',       ''                      };
        header_comment  ={'',       ''                      };

        txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
        txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
        txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';

        dlmwrite(filename,txt1,'delimiter','');
        dlmwrite(filename,txt2,'delimiter','','-append');
        dlmwrite(filename,txt3,'delimiter','','-append');
        dlmwrite(filename,[dz' p_dz'],'delimiter','\t','-append');

    % statistic 3 : curvature distribution absolute + sign
        filename = fullfile(dest,sprintf('%s%s',prefix,'statistic3.dat'));

        header_longname ={'mean curvature (absolute value)',    'total curvature (absolute value)',   'curvature sign',     'curvature distribution'};
        header_unit     ={sprintf('1/%s',unit),                  sprintf('1/%s',unit),                 '',                          unit};
        header_comment  ={'',                                    '',                                  '',                          ''};

        txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
        txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
        txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';

        dlmwrite(filename,txt1,'delimiter','');
        dlmwrite(filename,txt2,'delimiter','','-append');
        dlmwrite(filename,txt3,'delimiter','','-append');
        dlmwrite(filename,[abs(d2z') 2*abs(d2z') sign(d2z') p_d2z'],'delimiter','\t','-append');

    % Export film position data to "film_position.dat" file

    filename = fullfile(dest,sprintf('%s%s',prefix,'film_position.dat'));
    
    longname_text   = cell(length(data1.theta_vec),1); longname_text(:) = {''}; 
    unit_text       = cell(length(data1.theta_vec),1); unit_text(:)     = {unit};
    comment_text    = strread(num2str(atand(data1.theta_vec)),'%s');
    
    header_longname ={'total curvature'           longname_text{:}};
    header_unit     ={sprintf('1/%s',unit)        unit_text{:}};
    header_comment  ={''                          comment_text{:}};
    
    txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
    txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
    txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';
    
    dlmwrite(filename,txt1,'delimiter','');
    dlmwrite(filename,txt2,'delimiter','','-append');
    dlmwrite(filename,txt3,'delimiter','','-append');
    dlmwrite(filename,[2*data1.H_vec' data1.AVG_FILM_POSITION'],'delimiter','\t','-append');

    % Export volume data to "volume.dat" file

    filename = fullfile(dest,sprintf('%s%s',prefix,'volume.dat'));
    
    longname_text   = cell(length(data1.theta_vec),1); longname_text(:) = {''}; 
    unit_text       = cell(length(data1.theta_vec),1); unit_text(:)     = {sprintf('%s^3',unit)};
    comment_text    = strread(num2str(atand(data1.theta_vec)),'%s');
    
    header_longname ={'tota curvature'           longname_text{:}};
    header_unit     ={sprintf('1/%s',unit)       unit_text{:}};
    header_comment  ={''                         comment_text{:}};
    
    txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
    txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
    txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';
    
    dlmwrite(filename,txt1,'delimiter','');
    dlmwrite(filename,txt2,'delimiter','','-append');
    dlmwrite(filename,txt3,'delimiter','','-append');
    dlmwrite(filename,[2*data1.H_vec' data1.VOLUME'],'delimiter','\t','-append');

    % write phase diagram data to "pd.dat" file                
    
    filename = fullfile(dest,sprintf('%s%s',prefix,'phase_diagram.dat'));
    
    header_longname ={'theta [rad]',    'theta [deg]',      'total curvature',            'interface_position',   'jump_height'   };
    header_unit     ={'rad',            'deg',              sprintf('1/%s',unit),   unit,                   unit            };
    header_comment  ={'',               '',                 '',                     '',                     ''              };

    txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
    txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
    txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';
    
    dlmwrite(filename,txt1,'delimiter','');
    dlmwrite(filename,txt2,'delimiter','','-append');
    dlmwrite(filename,txt3,'delimiter','','-append');
    
    dlmwrite(filename,[phase_diagram.theta' atand(phase_diagram.theta)' 2*phase_diagram.curvature phase_diagram.position phase_diagram.amplitude],'delimiter','\t','-append');

    % write percolation diagram data to "percolation.dat" file                
    
    filename = fullfile(dest,sprintf('%s%s',prefix,'percolation_diagram.dat'));
    
    header_longname ={'theta [rad]',    'theta [deg]',      'critical curvature'};
    header_unit     ={'rad',            'deg',              sprintf('1/%s',unit)};
    header_comment  ={'',               '',                 ''};

    txt1=sprintf('%s\t',header_longname{:}); txt1(end)='';
    txt2=sprintf('%s\t',header_unit{:}); txt2(end)='';
    txt3=sprintf('%s\t',header_comment{:}); txt3(end)='';
    
    dlmwrite(filename,txt1,'delimiter','');
    dlmwrite(filename,txt2,'delimiter','','-append');
    dlmwrite(filename,txt3,'delimiter','','-append');
    
    dlmwrite(filename,[data1.theta_vec' atand(data1.theta_vec)' 2*data1.Hp'],'delimiter','\t','-append');

    % Export example of graphical resolution to "example.dat" file

    %to do
    %xlswrite('example.xls', [z' smooth(Lambda,0.1,'lowess') smooth(2*data1.W*curvature,0.1,lowess)],'grapres example');