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

% Load some packages
pkg load signal

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

% start loops
flow = loops.Flow;
pressure = loops.Pressure;

% time for plotting
time = (1:size(flow))*(1/Hz);

clc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Find the phase difference between pressure and
% flow by splitting into half breaths.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filter the pressure a little, because it's MESSY
pressure = sgolayfilt(pressure);

% Find the indices of insp start and exp start
flow_splits = splitBreaths(flow);
pressure_splits = splitBreaths(-pressure);

% Start and end stored in matrix
START = 1;
END = 2;

% plot pressure and reduced flow with split points
figure(1)
hold on
plot(-pressure, '.-b')
plot(flow/10, '.-k')
for value = 1:length(pressure_splits)
    plot(pressure_splits(START, value), -pressure(pressure_splits(START, value)), '.g')
    plot(pressure_splits(END, value), -pressure(pressure_splits(END, value)), '.r')
end
for value = 1:length(flow_splits)
    plot(flow_splits(START, value), flow(flow_splits(START, value))/10, '.c')
    plot(flow_splits(END, value), flow(flow_splits(END, value))/10, '.m')
end
legend('Pressure', 'Flow x 0.1')
grid minor
hold off

close(1)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% First method:
% Look at the offset in data points
%
% Here we are only looking at small/normal
% breaths. This keeps the acceleration of
% thew air at a low rate so the compliance
% can be found.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Start of inh. and exh. stored in matrix
INHALE = 1;
EXHALE = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ONLY LOOKING AT INHALATION AT THE MOMENT %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get each pair from pressure and flow
pressure_pair = 1
% start and end indices of inspiration
pressure_start = pressure_splits(INHALE, pressure_pair);
pressure_end = pressure_splits(EXHALE, pressure_pair);

flow_pair = 1
% start and end indices of inspiration
flow_start = flow_splits(INHALE, flow_pair);
flow_end = flow_splits(EXHALE, flow_pair);

while(pressure_pair < length(pressure_splits) && flow_pair < length(flow_splits))
    % Work out length of whole breath assuming
    % inh ans exh are the same length and shape
    % for a breath
    breath_length = 2*(pressure_end - pressure_start);
    % Check that the flow and pressure begin within
    % a half breath of eachother (with compliance it
    % can't be further than a half breath out)
    [pressure_pair, flow_pair] = checkIndicesAreForSameBreath(pressure_pair, pressure_splits, ...
                                flow_pair, flow_splits, breath_length)
    % The function will return non-indexed values
    % if no pair was found.
    if(pressure_pair > length(pressure_splits) || flow_pair > length(flow_splits))
        disp('broke');
        break;
    end
    pressure_pair = pressure_pair + 1;
    flow_pair = flow_pair + 1;
end
% throw away lower pair and pick up next pair
% Work out datapoint difference between pair start and end
% Work out radian per data point -> radian offset
% Get frequency (1/T)
% Find the resistance for spirometer pressure dropper
% The next bit is haaaaard =/


