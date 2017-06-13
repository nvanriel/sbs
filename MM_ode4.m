function [f] = MM_ode4(t,x,par,inputfile);
% FUNCTION MM_ode4 Full kinetic model of irreversible enzyme reaction with 
% inflow and outflow. Variable input.
%
%   [f] = MM_ode4(t,x)
%
%  Inputs: t - time (seconds)
%          x - vector of concentrations [a,b,c] (M)
%          par - vector of parameters [kp1, km1, kp2, E0, kp3]
%          input - function name in which input is defined
%  Outputs: f - vector of time derivatives
%              [da/dt,db/dt,dc/dt]      

%State variables
a = x(1);   %substrate
b = x(2);   %product
c = x(3);   %enzyme complex

%Parameters
kp1 = par(1);   %kp1 - forward rate constant (M^{-1} sec^{-1})
km1 = par(2);   %km1 - reverse rate constant (sec^{-1})
kp2 = par(3);   %kp2 - forward rate constant (sec^{-1})
E0  = par(4);   %E0 - total enzyme concentration (M)
kp3 = par(5);   %kp3 - rate constant product outflux (sec^{-1})

%Input
u = inputfile(t);

%ODE's
f(1) = u + (km1 + kp1*a)*c - kp1*E0*a;
f(2) = kp2*c - kp3*b;
f(3) = kp1*a*(E0 - c) - (km1 + kp2)*c;

f(4) = u;

f = f(:);   %Output needs to be a column vector
