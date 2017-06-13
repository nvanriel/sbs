function MM_GUI2
%Example of GUI with multiple sliders to change parameters in kinetic model 
%of irreversible enzymatic reaction
close all

% Default parameter values:
kp1 = 1000; %kp1 - forward rate constant (M^{-1} sec^{-1})
km1 = 1.0;  %km1 - reverse rate constant (sec^{-1})
kp2 = 0.1;  %kp2 - forward rate constant (sec^{-1})
E0  = 1e-4; %E0 - total enzyme concentration (M)
global par
par = [kp1, km1, kp2, E0];

%Create a figure for the GUI and configure the axes for displaying the simulation.
f = figure('Position',[400 50 560 600]); %[left, bottom, width, height] 

ax = axes('Parent',f,'position',[0.13 0.55 0.77 0.4]);    %[left, bottom, width, height] 
%h = plot(ax,[]);
set(gca,'XLim',[0 1000],'YLim',[0,1.2]);
xlabel('Time (s)'); ylabel('(mM)')
title('parameter sensitivity')
hold on
vertposition = 200; %vertical position of first slider
maxfold = 5;    %max fold increase of parameter value

%% slider 1 - Add the slider and slider label text to the figure.
b = uicontrol('Parent',f,...
    'Style','slider',...
    'Position',[81,vertposition,419,23],...   %
    'value',kp1, 'min',0, 'max',maxfold*kp1,...
    'Callback', @slider1_callback );
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,vertposition,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,vertposition,50,23],...
                'String',num2str(maxfold*kp1),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,vertposition-25,100,23],...
                'String','kp1','BackgroundColor',bgcolor);

%% slider 2
vertposition = vertposition-50;
b = uicontrol('Parent',f,...
    'Style','slider',...
    'Position',[81,vertposition,419,23],...
    'value',km1, 'min',0, 'max',maxfold*km1,...
    'Callback', @slider2_callback );
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,vertposition,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,vertposition,50,23],...
                'String',num2str(maxfold*km1),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,vertposition-25,100,23],...
                'String','km1','BackgroundColor',bgcolor);
 
%% slider 3
vertposition = vertposition-50;
b = uicontrol('Parent',f,...
    'Style','slider',...
    'Position',[81,vertposition,419,23],...
    'value',kp2, 'min',0, 'max',maxfold*kp2,...
    'Callback', @slider3_callback );
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,vertposition,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,vertposition,50,23],...
                'String',num2str(maxfold*kp2),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,vertposition-25,100,23],...
                'String','kp2','BackgroundColor',bgcolor);

%% slider 4
vertposition = vertposition-50;
b = uicontrol('Parent',f,...
    'Style','slider',...
    'Position',[81,vertposition,419,23],...
    'value',E0, 'min',0, 'max',maxfold*E0,...
    'Callback', @slider4_callback );
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,vertposition,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,vertposition,50,23],...
                'String',num2str(maxfold*E0),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,vertposition-25,100,23],...
                'String','E0','BackgroundColor',bgcolor);

text(0,-0.30, 'Click on any slider to simulate with default parameters')
text(0,-0.37, 'Move sliders to change values of corresponding parameters')


function slider1_callback(source, callbackdata, par)
%source 
%callbackdata
global par
%par = [kp1, km1, kp2, E0];
par(1) = source.Value;
simul(par);

function slider2_callback(source, callbackdata, par)
global par
par(2) = source.Value;
simul(par);

function slider3_callback(source, callbackdata, par)
global par
par(3) = source.Value;
simul(par);

function slider4_callback(source, callbackdata, par)
global par
par(4) = source.Value;
simul(par);


function simul(par)
% par: Parameter values
% par = [kp1, km1, kp2, E0, kp3];

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
