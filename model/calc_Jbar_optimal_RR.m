function Jbar_optimal = calc_Jbar_optimal_RR(theta,p)
%JBAR_OPTIMAL computes optimal Jbar for a given parameter combination and
%probe probability, according to the resource rational model (van den Berg
%& Ma, 2018). 
%
% ================== INPUT VARIABLES ===========
% 
%   THETA: parameter vector
%         tau: spread parameter of gamma function
%      lambda: weight for neural cost
%        beta: exponent for behavioral cost
% 
%   P: probability of being probed
% 
% ================== OUPUT VARIABLE ==================
% 
%   JBAR_OPTIMAL: optimal Jbar for given 
%
% =====================
%      Aspen Yoo
%   aspen.yoo@nyu.edu
% ======================

gvar = loadvar('gvar');

tau=theta(1);
lambda=theta(2);
beta=theta(3);

Jbar_optimal = fminsearch(@(pars) cost_function(pars,tau,lambda,beta,p,gvar), 1);

% This function returns expected total cost for a given set of parameters and p_i value
function E_C_total=cost_function(Jbar,tau,lambda,beta,p_i,gvar)
J = discretize_gamma(Jbar,tau,gvar.n_gamma_bins);
kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap)));

VM_x = linspace(-pi,pi,gvar.n_VM_bins);
VM_x = VM_x(2:end)-diff(VM_x(1:2))/2;
VM_y = exp(bsxfun(@times,kappa',cos(VM_x)));
VM_y = bsxfun(@rdivide,VM_y,sum(VM_y,2));

C_behavioral = abs(VM_x).^beta;
E_C_behavioral = mean(sum(bsxfun(@times,VM_y,C_behavioral),2));
E_C_total = p_i*E_C_behavioral + lambda*Jbar;

% This function discretizes a gamma distribution and returns bin centers
function bins = discretize_gamma(Jbar,tau,nbins)
X = linspace(0,1,nbins+1);
X = X(2:end)-diff(X(1:2))/2;
warning off
bins = gaminv(X,Jbar/tau,tau);
warning on