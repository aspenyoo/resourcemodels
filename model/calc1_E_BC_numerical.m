function E_C_behavioral = calc1_E_BC_numerical(Jbar, tau, beta)
%CALC_E_BC calculates expected behavioral cost given parameters

[JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');

% get Jpdf: p(J|Jbar,tau). 500 x 1 vector
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
Jpdf = Jpdf./sum(Jpdf);

d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);

% get dpdf (marginalize over J)
dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
dpdf = dpdf./sum(dpdf);

C_behavioral = dVec.^beta;
E_C_behavioral = sum(C_behavioral .* dpdf); % expected cost (marginalized over actual error)
