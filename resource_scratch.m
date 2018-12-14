

%% for a given Jbar_total, how the proportion allocation differ

clear all;

nSteps = 10;
p2Vec = linspace(0,0.5,nSteps);

for istep = 1:nSteps
    exppriorityVec = [0.5 p2Vec(istep) 0.5-p2Vec(istep)];
    
    
    
end

% get Jbar_total from RR model