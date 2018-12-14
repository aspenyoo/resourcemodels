function varargout = loadvar(varargin)
varvec = nan(1,nargin);
for iarg = 1:nargin
    varvec(iarg) = ischar(varargin{iarg});
end
varvec = find(varvec);
% nvars = sum(isstr(varargin));
nvars = length(varvec);
varargout = cell(1,nvars);
for ivar = 1:nvars
    var = varargin{varvec(ivar)};
    switch var
        case 'JVec'
            Jbar = varargin{varvec(ivar)+1}{1};
            tau = varargin{varvec(ivar)+1}{2};
            nJSamp = 500;
            
%             JVec = linspace(1e-10,10*Jbar,nJSamp);
            X = [.001 0.99];
            warning off
            JVec = gaminv(X,Jbar/tau,tau);
            warning on
            JVec = linspace(JVec(1),JVec(2),nJSamp)';
            
            varargout{ivar} = JVec;
        case 'rVec'
%             nRs = varargin{2*ivar};
            % radius stuff
            nRs = 1000;
            rVec = linspace(0,500,nRs);
            varargout{ivar} = rVec;
            
        case 'gvar'
            % a structure for easy passing to functions for RR model
            gvar.kmap = [linspace(0,10,250) linspace(10.001,20000,500)];
            gvar.Jmap = gvar.kmap.*besseli(1,gvar.kmap,1)./besseli(0,gvar.kmap,1); % mapping from kappa to J (see Appendix 1)
            gvar.n_gamma_bins = 500;        % number of bins to use to discretize Gamma distribution over J when computing model predictions
            gvar.n_VM_bins = 360;          % number of bins to use to discretize Von Mises distribution over estimation error when computing model predictions
            gvar.n_initpars_try = 25;      % number of random parameter vectors to try for initialization of optimizer
            %             gvar.uPi = unique(data.p_i)';  % unique values of p_i in this dataset
            varargout{ivar} = gvar;

    end
end