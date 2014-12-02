function [] = plotSolution(sigma1,sigma2,L,W,z,angle,curvature)
% Plot graphical resolution of desorption isotherm
%
% Args:
%   - sigma1, sigma2, L, W: functions defining the main equation of
%                           Wenzel Prewetting model (see main script)
%   - z: z coordinates (on which sigma1/2, L and W are indexed)
%   - angle: contact angle in degree
%   - curvature: total curvature
%
% Returns:
%   - a plot of the graphical resolution

theta  = tand(angle);   % angle in radian
Lambda = (sigma2./sigma1-theta).*L;

figure; hold on; subplot(1,2,1)
plot(z,smooth(Lambda,0.1,'lowess'),'ob',z,Lambda,'-b',z,2*W*curvature,'-r','erasemode','background');
subplot(1,2,2);
plot(z,smooth(Lambda,0.1,'lowess')'-2.*W*curvature,'-or','erasemode','background');
hold off;

end

