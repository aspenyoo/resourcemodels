%% gamma distribution

clear all
Jbar = .2; 
tau = .2; 

JVec = loadvar('JVec',{Jbar,tau});
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value

plot(JVec,Jpdf);

%% for a given Jbar_total, how the proportion allocation differ

clear all;

nSteps = 10;
p2Vec = linspace(0,0.5,nSteps);

for istep = 1:nSteps
    exppriorityVec = [0.5 p2Vec(istep) 0.5-p2Vec(istep)];
    
    
    
end

% get Jbar_total from RR model

%% comparing analytical and numberical calculations of behavioral loss

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