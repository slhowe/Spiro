% Author: Sarah Howe
% Written in octave 4.01
% May not be compatible with matlab

% Clean up
clc
close all
clear

% Automatic debugging
debug_on_interrupt(0);
debug_on_warning(1);
debug_on_error(1);

% Load some data
load SpirometryData.mat;
dataset = 1;
loops = data(dataset).Loops;
band = data(dataset).Banding;
normal = data(dataset).Normal;
inflated = data(dataset).Inflated;

% sampling frequency 125 Hz
Hz = 125;

% Atmospheric pressure in kpa
p_atm = 101.35; %kPa

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Treating the data as from an EoR-circuit
%
% Flow is current
% Pressure is voltage
% Resistance is resistance
% Compliance is capacitance
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% start loops
flow = loops.Flow;
pressure = loops.Pressure;

% time for plotting
time = (1:size(flow))*(1/Hz);

clc

figure(1)
hold on
plot(0.1*flow(1:2000), 'k')
plot(-pressure(1:2000), 'm')
grid minor
legend("0.1 x flow", "Pressure")
hold off

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Find the phase difference between pressure and
% flow. Capacitance (lung compliance) will cause
% a phase shift proportional to its size
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
i = splitBreaths(pressure)
