function E_C_behavioral = calc1_E_BC_numerical_rvdb(Theta)
%CALC_E_BC calculates expected behavioral cost for one priorities'
%parameters, using Ronald van den Berg's way of numerically integrating

% getting parameters
Jbar = Theta(1);
tau = Theta(2);
% lambda = Theta(end-1);
gamma = Theta(end);

% loading variables
[dVec,gvar] = loadvar('rVec','gvar');
JVec = discretize_gamma(Jbar,tau,gvar.n_gamma_bins); % particular J values based on Jbar and tau

d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);

C_behavioral = dVec.^gamma;
E_C_behavioral = mean(sum(bsxfun(@times,d_given_J,C_behavioral),2)); % expected cost (marginalized over actual error)

% This function discretizes a gamma distribution and returns bin centers
function bins = discretize_gamma(Jbar,tau,nbins)
X = linspace(0,1,nbins+1);
X = X(2:end)-diff(X(1:2))/2;
warning off
bins = gaminv(X,Jbar/tau,tau)';
warning on