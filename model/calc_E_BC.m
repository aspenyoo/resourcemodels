function E_C_total = calc_E_BC(Theta,allocatedpriorityVec,exppriorityVec)
%CALC_E_BC calculates expected behavioral cost given parameters



% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
beta = Theta(end);

nPriorities = length(allocatedpriorityVec);

E_C_total = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);
    p_i = exppriorityVec(ipriority);

    E_C_total=E_C_total + p_i*calc1_E_BC([Jbar Theta(2:end)]);
end