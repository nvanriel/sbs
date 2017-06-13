function MM_GUI
%Example of GUI with a slider to change a parameter in kinetic model of 
%irreversible enzymatic reaction
close all

% Default parameter values:
kp1 = 1000; %kp1 - forward rate constant (M^{-1} sec^{-1})
km1 = 1.0;  %km1 - reverse rate constant (sec^{-1})
kp2 = 0.1;  %kp2 - forward rate constant (sec^{-1})
E0  = 1e-4; %E0 - total enzyme concentration (M)
kp3 = 0.01;  %kp3 - rate constant product outflux (sec^{-1})


%Create a figure for the GUI and configure the axes for displaying the simulation.
f = figure;
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
%h = plot(ax,[]);
set(gca,'XLim',[0 1000],'YLim',[0,1.2]);
xlabel('Time (s)'); ylabel('(mM)')
%title('irreversible enzymatic reaction')
hold on

%Add the slider and slider label text to the figure.
b = uicontrol('Parent',f,...
    'Style','slider',...
    'Position',[81,54,419,23],...
    'value',kp1, 'min',0, 'max',1e4,...
    'Callback', @simul);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','1e4','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','kp1','BackgroundColor',bgcolor);

function simul(source, callbackdata)
%source 
%callbackdata

% Parameter values:
%kp1 = 1000; %kp1 - forward rate constant (M^{-1} sec^{-1})
kp1 = source.Value;
km1 = 1.0;  %km1 - reverse rate constant (sec^{-1})
kp2 = 0.1;  %kp2 - forward rate constant (sec^{-1})
E0  = 1e-4; %E0 - total enzyme concentration (M)
kp3 = 0.01;  %kp3 - rate constant product outflux (sec^{-1})

par = [kp1, km1, kp2, E0, kp3];

% Initial Conditions:
x0 = [0.001 0 0];

% Simulation settings:
tspan = [0 5000];    %(s)
odeoptions = [];    %use defaults
[t,x] = ode15s(@MM_ode,tspan,x0,odeoptions, par);

% Plot results:
plot(t,x(:,1)*1e3,'b', t,x(:,2)*1e3,'r', t,x(:,3)*1e3,'y');
legend('a','b','c')

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
