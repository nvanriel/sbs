%MM_script5 Full kinetic model of irreversible enzyme reaction with 
%inflow and outflow.
%Simulation with variable step integration method (ode15s).
%Parameters and input are defined outside the ODE file

% Parameter values:
kp1 = 1000; %kp1 - forward rate constant (M^{-1} sec^{-1})
km1 = 1.0;  %km1 - reverse rate constant (sec^{-1})
kp2 = 0.1;  %kp2 - forward rate constant (sec^{-1})
E0  = 1e-4; %E0 - total enzyme concentration (M)
kp3 = 0.01;  %kp3 - rate constant product outflux (sec^{-1})
par = [kp1, km1, kp2, E0, kp3];
%Input
inputfile = @MM_data;
% Initial Conditions:
x0 = [0.001 0 0 0];
% Integrate ODE:
tspan = [0 5000];    %(s)
odeoptions = [];    %use defaults
[t,x] = ode15s(@MM_ode4,tspan,x0,odeoptions, par,inputfile);
% Plot results:
figure; plot(t,x(:,1)*1e3,t,x(:,2)*1e3,t,x(:,3)*1e3);
xlabel('Time (s)'); ylabel('(mM)')
legend('a','b','c')
title('Model with measured influx as input')

figure; plot(t(2:end),diff(x(:,4))./diff(t))
