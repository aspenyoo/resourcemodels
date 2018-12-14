function E_BC = calc1_E_BC_analytical(Jbar, tau, beta)
%CALC_E_BC_ANALYTICAL calculates the expected behavioral cost given JBAR,
%TAU, and BETA. 
% 
% aspen, write conditions of the analytical solution

k = Jbar/tau;
E_BC = exp(gammaln(beta/2 + 1) + gammaln(k-(beta/2)) - gammaln(k)) .* (2/tau)^(beta/2);



