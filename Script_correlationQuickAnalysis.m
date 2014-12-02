% Compute radial autocorrelation
[R, rho] = radial_ACF(X,Y,Z,400);  % Compute radial ACF

% Get Root Mean Square roughness
sigma = sqrt(R(1));

% Gaussian fit (to estimate correlation length)
idx1 = ~isnan(rho) & ~isnan(R);
f = fit(rho(idx1)',R(idx1)'/sigma^2,'exp(-(x/T)^2)');
coefs1 = coeffvalues(f); cl = coefs1(1);

% polynomial fit (to get Dx coefficients)
idx2 = ~isnan(rho) & ~isnan(R) & rho<cl;
polystring = '1 - D2*x^2 - D4*x^4';
poly = fit(rho(idx2)',R(idx2)',polystring);
coefs2 = coeffvalues(poly); D2 = coefs2(1); D4 = coefs2(2);

%plot
figure; hold on;
plot(rho,R/sigma^2,'-');
set(gca,'XScale','log');
plot(f)
plot(cl,exp(-1),'or')