function PC_circles = Pair_correlation_ising(initials, domain)

% For fitting Pair correlation functions generated by get_autocorr.m
%
% domain = xdata used in original data

% Variables used:
% sig_s = PSF standard deviation
% rho_p = mean concentration of protein

% A_ising = pre-exponential amplitude 
% xi_ising = correlation length


sig_s = initials(1);
rho_p = initials(2);
A_ising = initials(3);
xi_ising = initials(4);


% For convolution purposes, spread out domain in both directions by a good
% amount.

min_dom = min(domain);
max_dom = max(domain);
dom_step = domain(2) - domain(1);

%r = (-5*max_dom):(dom_step):(5*max_dom);

r = -5*max_dom:(dom_step):(5*max_dom);

%%%% Elementary functions

g_stoch = (1/(4*pi*rho_p*(sig_s^2))).*exp(-(r.^2)/(4*sig_s^2));

g_PSF = (1/(4*pi*(sig_s^2))).*exp(-(r.^2)/(4*pi*sig_s^2));

g_ising = (A_ising*exp(-(r/xi_ising)).*abs((r.^(-1/4)))) + 1;

% This explodes at 0.  fix this.
% Set everything on the negative side to 1.
g_ising(le(r, 0)) = 1;

g_meas_ising= 1 + conv((g_ising-1), g_PSF, 'same') + g_stoch; 

    
%%% Cut assmbled functions down to size of original domain

r_start = find(min(abs(r - min_dom)) == abs((r - min_dom)));

r_end = find(min(abs(r - max_dom)) == abs((r - max_dom)));
%r_plot = r(r_start:r_end);

PC_circles = g_meas_ising(r_start:r_end);


