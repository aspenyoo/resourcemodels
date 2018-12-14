function E_C_behavioral = calc1_E_BC_halfnumerical_ahy(Theta)
%CALC_E_BC calculates expected behavioral cost given parameters

% getting parameters
Jbar = Theta(1);
tau = Theta(2);
beta = Theta(end);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

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

% d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);

% get dpdf (marginalize over J)
% dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
% dpdf = dpdf./sum(dpdf);

% C_behavioral = dVec.^gamma;
% E_C_behavioral = sum(C_behavioral .* dpdf); % expected cost (marginalized over actual error)
