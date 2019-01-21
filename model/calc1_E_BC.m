function E_BC = calc1_E_BC(Jbar, tau, beta)
%CALC_E_BC calculates expected behavioral cost given parameters

if (Jbar == 0); E_BC = Inf; return, end

k = Jbar/tau;
if (beta > -2) && ((k-beta/2) > 0) % if both analytical assumptions are satisfied
%     sprintf('analytical')
    E_BC = calc1_E_BC_analytical(Jbar, tau, beta);
elseif (beta > -2) % if only the first part is satisfied, can numerically integrate over just J
%     sprintf('half numerical')
    E_BC = calc1_E_BC_halfnumerical(Jbar, tau, beta);
else % if neither are satisfied, integrate over both d and J numerically
%     sprintf('numerical')
    E_BC = calc1_E_BC_numerical(Jbar, tau, beta);
end
