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