function expectederror = calc_E_BC_numerical_ahy(Theta,allocatedpriorityVec,exppriorityVec)
%CALC_E_BC calculates expected behavioral cost given parameters



% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
lambda = Theta(end-1);
gamma = Theta(end);

gvar = loadvar('gvar');

nPriorities = length(allocatedpriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);
    p_i = exppriorityVec(ipriority);
    
    E_C_total=cost_function(Jbar,tau,lambda,gamma,p_i,gvar)';
    
    expectederror = expectederror + E_C_total;
    
%     % get x-axis of J and distance d 
%     [JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');
%     JVec = JVec';
%     
%     % get Jpdf: p(J|Jbar,tau). 500 x 1 vector
%     Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
%     Jpdf = Jpdf./sum(Jpdf);
%     
%     % get dpdf given J: p(d|1/J). will be nJs (500) x nds (500) vector
%     d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);
%     
%     % get dpdf (marginalize over J)
%     dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
%     dpdf = dpdf./sum(dpdf);
%     
%     % \int d^blah p(d) dd
%     expectederror = expectederror + priorityVec(ipriority).*sum((dVec.^gamma) .* dpdf);
    
end

% This function returns expected total cost for a given set of parameters and p_i value
function E_C_total=cost_function(Jbar,tau,lambda,gamma,p_i,gvar)
% p(J|\bar{J},\tau)
% JVec = discretize_gamma(Jbar,tau,gvar.n_gamma_bins); % particular J values based on Jbar and tau

[JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');

% get Jpdf: p(J|Jbar,tau). 500 x 1 vector
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
Jpdf = Jpdf./sum(Jpdf);

d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);

% get dpdf (marginalize over J)
dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
dpdf = dpdf./sum(dpdf);

C_behavioral = dVec.^gamma;
E_C_behavioral = sum(C_behavioral .* dpdf); % expected cost (marginalized over actual error)
E_C_total = p_i*E_C_behavioral + lambda*Jbar; % expected total cost

%--------

% kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap))); % corresponding kappas, using predefined Jbars and kappas in gvar
% 
% VM_x = linspace(-pi,pi,gvar.n_VM_bins); % equally spaced bins
% VM_x = VM_x(2:end)-diff(VM_x(1:2))/2;
% % p(error|J)
% VM_y = exp(bsxfun(@times,kappa',cos(VM_x))); % get von mises pdf of this
% VM_y = bsxfun(@rdivide,VM_y,sum(VM_y,2)); % normalized over samples
% 
% C_behavioral = abs(VM_x).^gamma; % behavioral cost for every amnt of error
% % \int error^\beta p(error|J) p(J|\bar{J},\tau) dJ d(error)
% % somehow by sampling J by inverse gamma, the marginalization over J is the
% % same as jss a mean.
% E_C_behavioral = mean(sum(bsxfun(@times,VM_y,C_behavioral),2)); % expected cost (marginalized over actual error)
% E_C_total = p_i*E_C_behavioral + lambda*Jbar; % epected total cost

% % This function discretizes a gamma distribution and returns bin centers
% function bins = discretize_gamma(Jbar,tau,nbins)
% X = linspace(0,1,nbins+1);
% X = X(2:end)-diff(X(1:2))/2;
% warning off
% bins = gaminv(X,Jbar/tau,tau)';
% warning on