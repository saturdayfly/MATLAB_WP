theta   = 16;                        % Contact angle value in degree
h0      = 1;                     % Film thickness at coexistence
hh      = 0;                    % Minimum allowed film thickness

% Compute parameters W0 and l0 for the short range potential function

dW      = 1 - cosd(theta);          % Work of adhesion
l0      = h0*2/(1+sqrt(5));
c       = (1+sqrt(5))/(1-sqrt(5));  % Constant
W0      = -c*exp(1/2+sqrt(5)/2)*dW; % Minimum potential value = W(h0)

clear h;
h = NaN(size(pressure));


for i = 1:1:length(pressure) % Compute film thickness for each value of pressure

    % Define short range potential function and its derivative

    %W  = @(x)   W0*(l0./x-1).*exp(-x./l0);
    dW = @(x)   W0.*exp(-x./l0).*(1/l0 - 1./x - l0./x^2) + pressure(i);

    h(i) = fzero(dW,h0); % Solve
    if h(i)<hh h(i)=hh; end
    
end

plot(pressure,h);
set(gca,'XScale','log');