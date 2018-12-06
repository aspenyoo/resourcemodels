function varargout = loadvar(varargin)
nvars = ceil(length(varargin)/2);
varargout = cell(1,nvars);
for ivar = 1:nvars
    var = varargin{2*ivar-1};
    switch var
        case 'JVec'
            Jbar = varargin{2*ivar}{1};
            tau = varargin{2*ivar}{2};
            nJSamp = 500;
            
            JVec = linspace(1e-10,10*Jbar,nJSamp);
            
%             nJSamp = 100;
%             xmin = Jbar;
%             xmax = Jbar;
% 
%                 % lower bound
%                 while gampdf(xmin,Jbar/tau,tau) > 1e-4
%                     if xmin <= 1
%                         xmin = 1e-10;
%                         break
%                     else
%                         xmin = xmin - 1;
%                     end
%                 end
%                 
%                 increment = 1;
%                 % upper bound
%                 while gampdf(xmax,Jbar/tau,tau) > 1e-4
%                     xmax = xmax + increment;
%                 end
%                 % decrease by smaller increments to make it more fine
%                 % grained
%                 increment = 0.1*increment;
%                 while gampdf(xmax,Jbar/tau,tau) < 1e-4
%                     xmax = xmax - increment;
%                     if xmax <= 0.1
%                         xmax = 0.05;
%                         break
%                     end
%                 end
% %             end
%             JVec = linspace(xmin,xmax,nJSamp);
% %             JVec = linspace(1e-5,20,nJSamp); % ASPEN: make sure range is reasonable for parameter range
            varargout{ivar} = JVec;
        case 'rVec'
%             nRs = varargin{2*ivar};
            % radius stuff
            nRs = 500;
            rVec = linspace(0,10,nRs);
            varargout{ivar} = rVec;
            
        case 'gvar'
            % a structure for easy passing to functions for RR model
            gvar.kmap = [linspace(0,10,250) linspace(10.001,20000,500)];
            gvar.Jmap = gvar.kmap.*besseli(1,gvar.kmap,1)./besseli(0,gvar.kmap,1); % mapping from kappa to J (see Appendix 1)
            gvar.n_gamma_bins = 50;        % number of bins to use to discretize Gamma distribution over J when computing model predictions
            gvar.n_VM_bins = 360;          % number of bins to use to discretize Von Mises distribution over estimation error when computing model predictions
            gvar.n_initpars_try = 25;      % number of random parameter vectors to try for initialization of optimizer
            %             gvar.uPi = unique(data.p_i)';  % unique values of p_i in this dataset
            varargout{ivar} = gvar;

    end
end