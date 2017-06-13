%MM_script2 Full kinetic model of irreversible enzyme reaction with 
%inflow and outflow..
%Simulation with variable step integration method (ode15s).

% Initial Conditions:
x0 = [0.001 0 0];   %[a,b,c] (M)
% Integrate ODE:
tspan = [0 1000];    %(s)
[t,x] = ode15s(@MM_ode2,tspan,x0);
% Plot results:
figure; plot(t,x(:,1)*1e3,t,x(:,2)*1e3,t,x(:,3)*1e3);
xlabel('Time (s)'); ylabel('(mM)')
legend('a','b','c')
title('Model with constant inflow')
