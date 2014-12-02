function output = computeIsotherms(p_f,sigma1,sigma2,L,W,z,H_max,theta_inc,nvalues_H,p,threshold,method,span,display)

% Compute desorption isotherms following the method described in :
%
% [1] S. Herminghauss, PRL, 109, 2012
%     "Universal Phase Diagram for Wetting on Mesoscale Roughness"
%
% Args: 
%   - sigma1     : vector of length N, average slope
%   - sigma2     : vector of length N, average square slope
%   - L          : vector of length N, length of contour line
%   - W          : vector of length N, area enclosed by the contour line
%   - z          : vector of length N, surface elevation
%   - dz         : vector of length M, surface slope
%   - H_max      : maximum mean curvature considered
%   - theta_inc  : increment in contact angle values
%   - nvalues_H  : number of curvature values considered
%   - p          : minimum value of mean curvature in power of 10.
%                  Example : p = -2 corresponds to Hmin = 1e-2
%   - threshold  : threshold z value used to compute percolation transition
%   - method     : method for smoothing (default 'lowess')
%   - span       : span for smoothing (default 0.1)
%   - display    : character string taking the value 'on' or 'off'
%                  enables to display graphical resolution while solving
%
% Returns :
%   - output : a structure containing the following fields :
%       .H_vec                 : Vector of curvature values used for the computation
%       .theta_vec             : Vector of theta values used for the computation
%       .AVG_FILM_POSITION     : Average position of the liquid film.
%                                   Matrix of size [length(H_vec) length(theta_vec)]
%       .VOLUME                : Adsorbed volume per unit area.
%                                   Matrix of size [length(H_vec) length(theta_vec)]
%       .thetaf                : critical contact angle
%       .Hp                    : mean curvature at which percolation occurs (connected to theta_vec)
%


%% Part 1 - Computation of adsorbed liquid

    % Compute maximum theta value theta_max = max(sigma2/sigma1) (see ref [1] eq. 19)
    
        theta_max =  max(smooth(z,sigma2./sigma1,span,method));   
        
    % Define contact angle and curvature vectors + allocate memory
            
        H_vec               = [0, logspace(p,log10(H_max),nvalues_H)];  % We use a logspaced vector here
        theta_vec           = tand(0:theta_inc:atand(theta_max));

        AVG_FILM_POSITION   = NaN(length(theta_vec),length(H_vec));
        VOLUME              = NaN(length(theta_vec),length(H_vec));
    
    % Go through all theta and H values
    
    w = waitbar(0,'Initializing waitbar...');
        
        for t = 1:1:length(theta_vec) % go through theta values
            
            waitbar(t/length(theta_vec),w,sprintf('Computing phase diagram. %.0f%% along...',t/length(theta_vec)*100))

            % Define Lambda + smooth (only changes when changing theta)
            Lambda  = smooth((sigma2./sigma1-theta_vec(t)).*L,span,method);
            
            for h = 1:1:nvalues_H % go through H values
                                
                % Define equation 12 of reference [1]
                data    = Lambda' - (2*W.*H_vec(h));
                
                %find where data crosses the 0 line with grad < 0 (gives the stable solution)
                idx = find(diff(data >= 0) & diff(data)<0, 1, 'last' ); 
                
                if ~isempty(idx)    % solution found --> use linear interpolation for better accuracy
                    
                    line_slope = (data(idx+1)-data(idx))./(z(idx+1)-z(idx));
                     AVG_FILM_POSITION(t,h) = z(idx)-data(idx)./line_slope;
                    
                    % Compute volume :
                    if idx > 1
                        VOLUME(t,h) = trapz(  z(1:max(idx)),( AVG_FILM_POSITION(t,h)-z(1:max(idx))).*p_f(1:max(idx)));
                    else % only 1 datapoint, cannot integrate, assume volume is zero
                        VOLUME(t,h) = 0;
                    end
                    
                else                % no solution found - surface is dry

                    AVG_FILM_POSITION(t,h) = min(z(:));
                    
                end

                % Plot resolution if display is on
                if strcmp(display,'on') % input argument display = 'on'
                    pause(0.5);
                    plot(z,Lambda,'ob',z,smooth(Lambda,span,'loess'),'-b',z,2*W*H_vec(h),'or','erasemode','background');
                    if ~isempty(idx)
                        hold on;
                        plot([ AVG_FILM_POSITION(t,h)  AVG_FILM_POSITION(t,h)], [0 max(Lambda)],'-g','LineWidth',2);
                        hold off;
                    end
                    drawnow;
                end
            end
        end
        
        close(w);

%% Part 2 - Compute percolation transition

    % Allocate
    Hp = NaN(size(theta_vec));  % Critical curvature at chich percolation occurs
       
    for t = 1:1:length(theta_vec) % go through theta

        temp = AVG_FILM_POSITION(t,:);
        [~,idx] = find(temp<=threshold,1,'first');
        
        if ~isempty(idx) % solution found
        Hp(t) = H_vec(idx);
        end       
        
    end    

%% Part 3 - Assign dat to output structure

    output.thetaf               = theta_max;
    output.H_vec                = H_vec;
    output.theta_vec            = theta_vec;
    output.AVG_FILM_POSITION    = AVG_FILM_POSITION;
    output.VOLUME               = VOLUME;
    output.Hp                   = Hp;

end