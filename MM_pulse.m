function u = MM_pulse(t)
%MM_pulse Pulse inflow for model of irreversible enzyme reaction.
%
%   u = MM_pulse(t)
%

if t<2000 | t>3000
    u = 1e-6;  %influx (M/sec)
else
    u = 5e-6;  %influx (M/sec)
end
