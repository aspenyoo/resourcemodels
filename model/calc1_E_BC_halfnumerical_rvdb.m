function E_C_behavioral = calc1_E_BC_halfnumerical_rvdb(theta)
%CALC_E_BC calculates expected behavioral cost for one priorities'
%parameters, using Ronald van den Berg's way of numerically integrating

% getting parameters
Jbar = theta(1);
tau = theta(2);
% lambda = Theta(end-1);
beta = theta(end);

% loading variables
[gvar] = loadvar('gvar');
JVec = discretize_gamma(Jbar,tau,gvar.n_gamma_bins); % particular J values based on Jbar and tau
k = Jbar/tau;

% assumption required for integration over d
gam2 = beta/2 + 1;
assert(gam2>0)

% get Jpdf: p(J|Jbar,tau). 500 x 1 vector
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
Jpdf = Jpdf./sum(Jpdf);

% integrating over some function of J, once analytically integrated over d
constant = 2^(beta/2)/tau^k*exp(gammaln(gam2)-gammaln(k));
E_C_behavioral = constant*mean((JVec).^(k-(beta/2)-1).*exp(-JVec./tau));

% 
% d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);
% 
% C_behavioral = dVec.^gamma;
% E_C_behavioral = mean(sum(bsxfun(@times,d_given_J,C_behavioral),2)); % expected cost (marginalized over actual error)

% This function discretizes a gamma distribution and returns bin centers
function bins = discretize_gamma(Jbar,tau,nbins)
X = linspace(0,1,nbins+1);
X = X(2:end)-diff(X(1:2))/2;
warning off
bins = gaminv(X,Jbar/tau,tau)';
warning on