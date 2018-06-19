% Pair autocorrelation for Gaussian peaks within a random field.
% Equation 2 of Sengupta, Nature Protocols, 2013

function gStoch = GaussianPairCorr(parameters, r)

sigma = parameters(1);
rho = parameters(2);

gStoch = (1./(4*pi*(sigma.^2).*rho)).*exp((-(r.^2))./(4.*sigma.^2)) + 1;