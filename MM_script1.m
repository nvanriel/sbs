%MM_script1 Full kinetic model of irreversible enzyme reaction.
%Simulation with Runge-Kutta integration method (ode45) and variable step 
%method (ode15s).

% Initial Conditions:
x0 = [0.001 0 0];   %(M)
% Integrate ODE:
tspan = [0 500];    %(s)
tic
[t,x] = ode45(@MM_ode1,tspan,x0);   %Runge-Kutta
toc
% Plot results:
figure; plot(t,x(:,1)*1e3,t,x(:,2)*1e3,t,x(:,3)*1e3);
xlabel('Time (s)'); ylabel('(mM)')
legend('a','b','c')
title('Closed model with ode45')

tic
[t,x] = ode15s(@MM_ode1,tspan,x0);  %variable step for stiff systems
toc
% Plot results:
figure; plot(t,x(:,1)*1e3,t,x(:,2)*1e3,t,x(:,3)*1e3);
xlabel('Time (s)'); ylabel('(mM)')
legend('a','b','c')
title('Closed model with ode15s')
