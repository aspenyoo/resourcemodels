% this script is a composite of random things involved in exploring,
% analyzing, and displaying data, model, and model fits. It is generally
% organized in terms of data related and model related sections, although
% each section is a bunch of random things I did. 

%% ==================================================================
%                       EXPERIMENT RELATED
% ===================================================================


%% testing out precue
clear all

priorityset = [0.5 0.25 0.125 0.125];

% open screen
screenNumber = max(Screen('Screens'));
[windowPtr] = Screen('OpenWindow', screenNumber, [0 0 0]);
Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[w, h]=Screen('WindowSize', screenNumber);  % Screen resolution
screenResolution = [w h];                 % Screen resolution
center = screenResolution/2;       % Screen center

% get stuff for aperature
apertureSize= 400;
apertureRect = [center-apertureSize center+apertureSize];

% draw stuffs
Screen('FillOval',windowPtr, 115, apertureRect);
% drawPrecue_old(windowPtr,priorityset)
drawPrecue(windowPtr,priorityset)

Screen('Flip',windowPtr);
pause;
sca
%% ==================================================================
%                       DATA RELATED
% ===================================================================

%% look at eye data

clear all; close all

subjid = 'AHY';
pricond = 1;

% put it in format
ifg_fn = '~/Documents/MATLAB/iEye_ts/examples/p_500hz.ifg';

ii_params = ii_loadparams; % load default set of analysis parameters, only change what we have to

ii_params.trial_end_value = 8;   % XDAT value for trial end
ii_params.drift_epoch = 1:5; % XDAT values for drift correction
ii_params.calibrate_epoch = 7 ;   % XDAT value for when we calibrate (feedback stim)
ii_params.calibrate_select_mode = 'last'; % how do we select fixation with which to calibrate?
ii_params.calibrate_mode = 'scale'; % scale: trial-by-trial, rescale each trial; 'run' - run-wise polynomial fit
ii_params.blink_window = [200 200]; % how long before/after blink (ms) to drop?
ii_params.plot_epoch = 1:8;  % what epochs do we plot for preprocessing?
ii_params.calibrate_limits = 2.5; % when amount of adj exceeds this, don't actually calibrate (trial-wise); ignore trial for polynomial fitting (run)
ii_params.ppd = 31.8578; % for scanner, 1280 x 1024 - convert pix to DVA

edf_prefix = sprintf('eyedata_%s_pricond%d_',subjid,pricond);

% files
root = sprintf('/Volumes/data/resourcemodels/output/%s/',subjid);
edf_files = dir(fullfile(root,sprintf('%s*.edf',edf_prefix)));
nFiles = length(edf_files);

% create empty cell array of all our trial data for combination later on
ii_trial = cell(nFiles,1);

for ifile = 1:nFiles
    ifile
    
    % what is the output filename?
    preproc_fn = sprintf('%s%s_preproc.mat',root,edf_files(ifile).name(1:(end-4)));
    
    
    [ii_data, ii_cfg, ii_sacc] = ii_preproc(fullfile(root,edf_files(ifile).name),ifg_fn,preproc_fn,ii_params);
    
%     if ifile == 1
%         % plot some features of the data
%         % (check out the docs for each of these; lots of options...)
%         ii_plottimeseries(ii_data,ii_cfg); % pltos the full timeseries
%         
%         ii_plotalltrials(ii_data,ii_cfg); % plots each trial individually
%         
%         ii_plotalltrials2d(ii_data,ii_cfg); % plots all trials, in 2d, overlaid on one another w/ fixations
%     end
    
    % score trials
    % default parameters should work fine - but see docs for other
    % arguments you can/should give when possible
    [ii_trial{ifile},~] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc);
    
end

ii_sess = ii_combineruns(ii_trial);
save(sprintf('%s%s_pricond%d_ii_sess.mat',root,subjid,pricond),'ii_sess','edf_files','ii_params')

%% look at error as a function of priority

clear all

subjid = 'AHY';
pricond = 1;

% load experimental design matrix and eye data
root = sprintf('/Volumes/data/resourcemodels/output/%s/',subjid);
load(sprintf('%s%s_pricond%d_ii_sess.mat',root,subjid,pricond),'ii_sess') % eye data
load(sprintf('%s%s_pricond%d_designMat.mat',root,subjid,pricond));%,'designMat','settings')

idx_targetpri = 5; % designMat column index of the target priority
priorityVec = unique(settings.prioritySets(pricond,:)); % current priority set
nPriorities = length(priorityVec);

% AHY (as of 2/29/2019) had some weird saving issues. these are adjustments
% so the data and designMat are matched
nTrialsPerRun = settings.nTrials/settings.nRuns;
switch pricond
    case 1
        % sessions 1-5: everything worked "correctly", but there was a bug 
        % in the experimental code where 12:12:XX trials were repeated
        dmat = designMat(1:12,:);
        for irun = 2:5
            dmat = [dmat; designMat(((irun-1)*nTrialsPerRun):(irun*nTrialsPerRun),:)];
        end
    case 2
        % session 3 is not included in here, bc it was only one trial,
        % which causes errors in preproc
        dmat = designMat(1:12,:);
        for irun = [2 4 5]
            dmat = [dmat; designMat(((irun-1)*nTrialsPerRun):(irun*nTrialsPerRun),:)];
        end
    case 3
        % session 3 combined session 2 and 3, so 24th trial was not
        % repeated, otherwise "correct" (has bug but otherwise correct)
        dmat = designMat(1:12,:);
        dmat = [dmat; designMat(((2-1)*nTrialsPerRun):(3*nTrialsPerRun),:)];
        for irun = 4:5
            dmat = [dmat; designMat(((irun-1)*nTrialsPerRun):(irun*nTrialsPerRun),:)];
        end
end
designMat = dmat;
targetpri = designMat(:,idx_targetpri);

data_bypri = cell(1,nPriorities);
for ipriority = 1:nPriorities
    priority = priorityVec(ipriority);
    
    idx_pri = targetpri == priority;
    data_bypri{ipriority} = ii_sess.f_sacc_err(idx_pri);
end

[M_error_bypri, sem_error_bypri] = deal(nan(1,nPriorities));
for ipriority = 1:nPriorities
    
    M_error_bypri(ipriority) = nanmean(data_bypri{ipriority});
    sem_error_bypri(ipriority) = nanstd(data_bypri{ipriority})./sqrt(sum(~isnan(data_bypri{ipriority})));
end

figure;
errorbar(M_error_bypri, sem_error_bypri)
set(gca,'XTick',1:nPriorities,'XTickLabel',priorityVec)
defaultplot

%%

% which_excl = [11 13 20 21]; % which indices to reject
% 11: drift correction too big
% 13: find fixations in this trial (something about channels, not done yet)
% 20: % no saccades identified in response epoch
% 21: % none of the identified saccades (>0) passed exclusion criteria for primary saccade
% 22: i_sacc error too high

%% look at error stuff

ii_sess.f_sacc_err
ii_sess.i_sacc_err
% 
% 
% priorityVec = [0.6 0.3 0.1];
% nPriorities = length(priorityVec);
% 
% [nTrials,pErr,fErr] = deal(nan(nSubj,nPriorities));
% 
%     % ==== GET EYE DATA =====
%     
%     % load filename
%     root = sprintf('/Volumes/data/Pri_quad/behavior/eyetracking/%s/',subjid);
%     load(sprintf('%s%s_ii_sess.mat',root,subjid))
%     
%     % trial exclusion
% %     which_excl = [11 13 20 21]; % which indices to reject
% %     use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
% use_trial = ~isnan(sum(ii_sess.f_sacc_err,2));
%     fprintf('keeping %0.03f%% of trials\n',mean(use_trial)*100);
%     
%     % get primary and final saccade error
%     primaryError = ii_sess.i_sacc_err(use_trial,:);
%     finalError = ii_sess.f_sacc_err(use_trial,:);
%     % target = ii_sess.targ(use_trial,:);
%     
%     % ====== PRIORITY DATA =========
%     
%     filename = sprintf('%s_outputmatrix.mat',subjid);
%     load(filename)
%     
%     % make sure that they are the same size
%     assert(size(use_trial,1) == size(outputMatrix,1))
%     
%     priorities = outputMatrix(use_trial,7);
%     
%     for ipriority = 1:nPriorities
%         priority = priorityVec(ipriority);
%         idx = priorities == priority;
%         
%         nTrials(isubj,ipriority) = sum(idx);
%         pErr(isubj,ipriority) = mean(primaryError(idx));
%         fErr(isubj,ipriority) = mean(finalError(idx));
%     end
%     
% end
% nTrials
% pErr
% fErr

%% ==================================================================
%                   MODEL RELATED
% ===================================================================

%% gamma distribution

clear all
Jbar = .2; 
tau = .2; 

JVec = loadvar('JVec',{Jbar,tau});
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value

plot(JVec,Jpdf);

%% comparing analytical and numberical calculations of behavioral loss
% 12/14/2018: code no longer useable

clear all

% theta = [3.5 0.5 0 1];
theta = [2 0.2 0 1];
exppriorityVec = [0.6 0.3 0.1];

nSteps = 10;
p2Vec = linspace(0,0.5,nSteps);

[analytical_loss, numerical_loss, numerical_spen, halfnumerical_ahy] = deal(nan(1,nSteps));
[a_time,b_time,c_time] = deal(0);
for istep = 1:nSteps
    allocatedpriorityVec = [0.5 p2Vec(istep) 0.5-p2Vec(istep)];
    
    tic; 
    analytical_loss(istep) = calc_expectederror_analytical([theta(1:2) theta(end)],allocatedpriorityVec,exppriorityVec);
    a_time = a_time+toc;
    
    tic;
    numerical_loss(istep) = calc_E_BC_numerical(theta, allocatedpriorityVec, exppriorityVec);
    b_time = b_time+toc;
    
    tic;
    numerical_spen(istep) = calc_E_BC_numerical_spen(theta, allocatedpriorityVec, exppriorityVec);
    c_time = c_time+toc;
end

figure; hold on
plot(analytical_loss,'r-')
plot([find(analytical_loss==min(analytical_loss)) find(analytical_loss==min(analytical_loss))],[0 10],'r-')

plot((numerical_loss)./16,'g--')
plot([find(numerical_loss==min(numerical_loss)) find(numerical_loss==min(numerical_loss))],[0 10],'g--')

plot(numerical_spen,'b:')
plot([find(numerical_spen==min(numerical_spen)) find(numerical_spen==min(numerical_spen))],[0 10],'b:')

%% plotting gamma functions with the estimates of behav loss for each
% 12/14/2018: code no longer useable

clear all

JbarVec = [0.01 0.05 0.1 0.2 0.5 1 2];
tau = 0.2;
beta = 1;

nJbars = length(JbarVec);

[analytical, halfnumerical_ahy, halfnumerical_rvdb,numerical_rvdb, numerical_ahy] = deal(nan(1,nJbars));
[a_time,hn_r_time,hn_a_time, n_r_time, n_a_time] = deal(0);

for ijbar = 1:nJbars
    Jbar = JbarVec(ijbar);
    theta = [Jbar tau beta];

    subplot(1,nJbars+1,ijbar)
    
    %aspen numerical
    JVec = loadvar('JVec',{Jbar,tau});
    Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    plot(JVec,Jpdf,'g.'); hold on
    
    % ronald numerical
    gvar = loadvar('gvar');
    JVec = discretize_gamma(Jbar,tau,gvar.n_gamma_bins); % particular J values based on Jbar and tau
    Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    plot(JVec,Jpdf,'b.');
    
    defaultplot


    tic; analytical(ijbar) = calc1_E_BC_analytical(theta); a_time = a_time + toc;
    tic; halfnumerical_ahy(ijbar) = calc1_E_BC_halfnumerical_ahy(theta); hn_a_time = hn_a_time + toc;
%     tic; halfnumerical_rvdb(ijbar) = calc1_E_BC_halfnumerical_rvdb(theta); hn_r_time = hn_r_time + toc;
%     tic; numerical_rvdb(ijbar) = calc1_E_BC_numerical_rvdb(theta); n_r_time = n_r_time + toc;
    tic; numerical_ahy(ijbar) = calc1_E_BC_numerical_ahy(theta); n_a_time = n_a_time + toc;
end

subplot(1,nJbars+1,nJbars+1); hold on
plot(analytical,'k-');hold on
% plot(halfnumerical_ahy,'r--')
% plot(halfnumerical_rvdb,'m--')
% plot(numerical_rvdb,'g--')
plot(numerical_ahy,'b:')
% plot(analytical./min(analytical),'k-');hold on
% plot(halfnumerical_ahy./min(halfnumerical_ahy),'r--')
% plot(halfnumerical_rvdb./min(halfnumerical_rvdb),'m--')
% plot(numerical_rvdb./min(numerical_rvdb),'g--')
% plot(numerical_ahy./min(numerical_ahy),'b:')

a_time,hn_r_time,hn_a_time, n_r_time, n_a_time

%% for a given parameter set, how Jbar_total differs with increasing p_i

clear all;

theta = [0.5 1 1];

nSteps = 15;
pVec = linspace(0,1,nSteps+1);
pVec = pVec(1:end-1) + diff(pVec(1:2));

optJbarVec = nan(1,nSteps);
for istep = 1:nSteps
    p_i = pVec(istep);
    
    optJbarVec(istep) = calc_Jbar_optimal_RR(theta, p_i);  
end

plot(pVec,optJbarVec,'k-')
xlabel('probe probability')
ylabel('resource allocated')
defaultplot


%% how optimal allocation changes with experimental probe probability
% for arbitrary parameter combination for Minimizing Error model

clear all
close all
clc

% theta
theta = [.4 0.075 0.5];

% ------ SET UP TERNARY PLOT -----
colorMat = [1 0 0; 0 0 1; 0 0 0];

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% get probe probabilities that will be investigated
x = roundn(0:0.01:1,-3);
nSteps = length(x);
[X,Y] = meshgrid(x,x);
X = X(:);
Y = Y(:);
idx = (X+Y) > 1;
X(idx) = [];
Y(idx) = [];
Z = 1-X-Y;

% get Jbar_totals for each p_i
optJbarVec = nan(1,nSteps);
for istep = 1:nSteps
    p_i = x(istep);
    
    optJbarVec(istep) = calc_Jbar_optimal_RR(theta, p_i);  
end

% substitute in everything
XX = roundn(X,-3); YY = roundn(Y,-3); ZZ = roundn(Z,-3);
for ix = 1:nSteps
    val = x(ix);
    
    XX(XX==val) = optJbarVec(ix);
    YY(YY==val) = optJbarVec(ix);
    ZZ(ZZ==val) = optJbarVec(ix);
end

JbartotalVec = sum([XX YY ZZ],2);

% % get indices of any X Y Z
% threshold = 0.01;
% idx = (X<=threshold);
% idx = idx | (Y<=threshold);
% idx = idx | (Z<=threshold);
% 
% X(idx) = [];
% Y(idx) = [];
% Z(idx) = [];
% JbartotalVec(idx) = [];

hold on
% get x and y of actual values
xx=0.5-X.*cos(pi/3)+Y./2;
yy=0.866-X.*sin(pi/3)-Y.*cot(pi/6)/2;
sz = 40;
colormap('parula')
s = scatter(xx,yy,sz,JbartotalVec,'filled');

%% for a given Jbar_total, how the proportion allocation differ
% in this section, I found that if I pick Jbar_total from the RR model 
% and use the ME model on that, it gives the same exact
% allocatedpriorityVec.

clear all;

nSteps = 10;
p1 = 0.5;
p2Vec = linspace(0.01,0.499,nSteps);
p3Vec = (1-p1) - p2Vec;
p2Vec, p3Vec

tau = 0.5; 
lambda = 2;
beta = 1;

% ================ RR model ================
% calculate optimal Jbar_total for p1
Jbar_p1 = calc_Jbar_optimal_RR([tau lambda beta],p1);

JbarMat = nan(nSteps,3);
JbarMat(:,1) = Jbar_p1;

% calculate Jbar_optimal for p2 and p3
for istep = 1:nSteps
    pVec = [p2Vec(istep) p3Vec(istep)];
    
    for ip = 1:2;
        pi = pVec(ip);
        
        JbarMat(istep,ip+1) = calc_Jbar_optimal_RR([tau lambda beta], pi);   
    end
end

% Jbar_total from RR model, for each exppriorityVec
Jbar_total = sum(JbarMat,2);

% pVec, for each condition
pVec_RR = bsxfun(@rdivide,JbarMat,Jbar_total);

% ========================== ME model ============================
% assuming same Jbar + exppriorityVec, get allocatedpriorityVec

pVec_ME = nan(nSteps,3);
for istep = 1:nSteps
    theta = [Jbar_total(istep) tau beta];
    exppriorityVec = [p1 p2Vec(istep) p3Vec(istep)];
    
    pVec_ME(istep,:) = calc_pVec_minerror(theta, exppriorityVec);
    
end

pVec_RR - pVec_ME

%% RR vs ME
clear all

% two experiment priorities that should result in very different 
nItems = 4;
exppriorityMat = [  0.6  0.3   0.1     0;
                    0.5 0.25 0.125 0.125;
                    0.5  0.3   0.2   0.0];
% exppriorityMat = [ones(1,nItems)./nItems;
%                    0.7 ones(1,nItems-1)*.3/(nItems-1);
%                    0.8 ones(1,nItems-1)*.2/(nItems-1);
%                    0.9 ones(1,nItems-1)*.1/(nItems-1);
%                    0.9 0.1 zeros(1,nItems-2)]; % with one probe probablity = 0, drastically different predictions from RR model
nExps = size(exppriorityMat,1); % number of "experiments"

% parameters chosen from realistic values in van den berg & ma, 2018.
tau = 0.4;
lambda = 0.075;
beta = 0.5; 

% calculate the Jbar_total for RR for each priorityVec
JbaroptMat = nan(nExps,nItems);
for ip = 1:numel(exppriorityMat);
    JbaroptMat(ip) = calc_Jbar_optimal_RR([tau lambda beta], exppriorityMat(ip));
end
JbarOptVec = sum(JbaroptMat,2);

pVec_RR = bsxfun(@rdivide, JbaroptMat, JbarOptVec);

% see how ME would predict across conditions, keeping Jbar_total fixed
pVec_ME = cell(1,nExps);
for iJbar = 1:nExps
    Jbar_total = JbarOptVec(iJbar); % assume Jbar_total is fixed
    theta = [Jbar_total tau beta];
    
    pVec_ME{iJbar} = nan(nExps,nItems);
    for iexp = 1:nExps
        exppriorityVec = exppriorityMat(iexp,:);
        pVec_ME{iJbar}(iexp,:) = calc_pVec_minerror(theta,exppriorityVec);
    end
end

pVec_RR
for iJbar = 1:nExps
    pVec_ME{iJbar}
end

%% calculate KL divergence

KLMat = nan(nExps);
for iJbar = 1:nExps
    KLMat(:,iJbar) = sum(pVec_RR.*log(pVec_RR) - pVec_RR.*log(pVec_ME{iJbar}),2);
end


% .02 for nItems = 4;