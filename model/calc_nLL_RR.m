function nLL = calc_nLL_RR(theta,data,exppriorityVec,fixparams)
% calc_nLL_RR(x,data,exppriorityVec,fixparams);
% 
%   ================= INPUT VARIABLES ================
%
%   THETA: parameter vector
%       parameter descriptions: 
%           TAU: second parameter of gamma noise distribution
%           LAMBDA: weight for neural cost
%           BETA: exponent for behavioral loss
% 
%   DATA: cell of length nPriorities. the ith element of DATA should be
%     all data corresponding to EXPPRIORITYVEC(i) condition. the first
%     column should contain the magnitude of errors and the second column,
%     if available, should contain the corresponding circle wager radius
%     size. 
% 
%   EXPPRIORITYVEC: vector of priority values (in decreasing order of
%     priority) e.g., [0.6 0.3 0.1]
% 
%   FIXPARAMS: (optional). 2 x (number of fixed parameters) matrix. fixed 
%     parameters, such that the first row corresponds to the index and 
%     second row corresponds to the value of the fixed parameter. 
%
%
%   ================= OUTPUT VARIABLES ================
% 
%   NLL: negative log-likelihood

% calc_nLL_RR calculates the negative log-likelihood of the resource
% rational model
% 
% Aspen Yoo
% aspen.yoo@nyu.edu
% 
% edited from script from van den Berg & Ma, 2018. 

tau=theta(1);
lambda=theta(2);
beta=theta(3);

gvar = loadvar('gvar');

nPs = length(exppriorityVec);
% p_resp = zeros(1,numel(data.p_i));
nLL = 0;
for ip = 1:nPs
    p = exppriorityVec(ip);
             
    % Compute optimal Jbar for this value of p_i
    Jbar_optimal = fminsearch(@(pars) cost_function(pars,tau,lambda,beta,p,gvar), 1);
    % Compute probability of the subject's estimation errors under this value of Jbar_optimal
    J = discretize_gamma(Jbar_optimal,tau,gvar.n_gamma_bins);
    kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap)));
    nLL = nLL - sum(log(max(mean(bsxfun(@times,1./(2*pi*besseli(0,kappa)),exp(bsxfun(@times,kappa,cos(data{ip})))),2),1e-3)));
%     p_resp(idx) = mean(bsxfun(@times,1./(2*pi*besseli(0,kappa)),exp(bsxfun(@times,kappa,cos(data.error(idx))))),2);
end
% nLL = -sum(log(max(p_resp,1e-3)));

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