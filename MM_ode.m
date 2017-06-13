function [f] = MM_ode(t,x,par)
% FUNCTION MM_ode1 Full kinetic model of irreversible enzymatic reaction.
%
%   [f] = MM_ode1(t,x)
%
%  Inputs:  t - time (seconds)
%           x - vector of concentrations (a,b,c) (M)
%  Outputs: f - vector of time derivatives
%               {da/dt,db/dt,dc/dt}      

%State variables
a = x(1);   %substrate
b = x(2);   %product
c = x(3);   %enzyme complex

%Parameters
kp1 = par(1);   %kp1 - forward rate constant (M^{-1} sec^{-1})
km1 = par(2);   %km1 - reverse rate constant (sec^{-1})
kp2 = par(3);   %kp2 - forward rate constant (sec^{-1})
E0  = par(4);   %E0 - total enzyme concentration (M)

%ODE's
f(1) = -kp1*E0*a + (km1 + kp1*a)*c;
f(2) = kp2*c;
f(3) = -(km1 + kp2)*c + kp1*a*(E0 - c);

f = f(:);   %Output needs to be a column vector
