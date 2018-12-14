function E_C_behavioral = calc1_E_BC_halfnumerical(Jbar, tau, beta)
%CALC_E_BC calculates expected behavioral cost given parameters


% assumption required for integration over d
gam2 = beta/2 + 1;
assert(gam2>0)

JVec = loadvar('JVec',{Jbar,tau});
deltax = diff(JVec(1:2));
k = Jbar/tau;

% get Jpdf: p(J|Jbar,tau). 500 x 1 vector
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
Jpdf = Jpdf./sum(Jpdf);

% integrating over some function of J, once analytically integrated over d
constant = 2^(beta/2)/tau^k*exp(gammaln(gam2)-gammaln(k));
E_C_behavioral = constant*sum((JVec).^(k-(beta/2)-1).*exp(-JVec./tau))*deltax;
