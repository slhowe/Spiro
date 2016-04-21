% Written by Sarah Howe in octave
% May not be compatible with matlab
%
% Playing with spirometry data

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
loops = data.Loops;
band = data.Banding;
normal = data.Normal;
inflated = data.Inflated;

% sampling frequency 125 Hz
Hz = 125;

% Atmospheric pressure in kpa
p_atm = 101.35; %kPa

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Treating the data as from an RC-circuit
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

%-----------------------------------------------
% Looks like an RC curve at end of forced exp?
% Let's fit a line to it :)
%-----------------------------------------------

% RC curve range
start = 1190;
stop = 2000;

figure(1)
hold on
plot(flow(start:stop), 'b', 'linewidth', 2)
xlabel("dataPoint")
ylabel("flow")
grid minor
hold off

% Least squares fit of curve
curveStart = start;
curveDataEnd = stop;

% set stopping point of fit at certain %age flow drop
% (things get over constrained below 3% and fit gets bad)
drop = 0.1*(flow(curveStart)-flow(curveDataEnd));

% find the index of the stopping point
index = 0;
stillLooking = 1;
for i = curveStart:curveDataEnd
    if(stillLooking)
        if(flow(i) > drop)
            index = i;
            stillLooking = 0;
        end
    end
end
if index == 0
    error("Percentage flow drop specified not found in range")
end
curveStop = index;

% set up matrices
measurements = log(-flow(curveStart:curveStop)); %flow flipped for nicer maths
one = ones(1, curveStop-curveStart+1);
times = -(time(curveStart:curveStop)-time(curveStart));

% OMG least squares!!!
results = [one', times']\measurements;

% extract info
startPoint = exp(results(1));
EoR = results(2)

% remake curve from info
times = -(time(curveStart:curveDataEnd)-time(curveStart));
newValues = -startPoint*exp(times*EoR); % flip flow to negative side

% compare new and old
figure(1)
hold on
plot(newValues, 'm')
hold off

%-----------------------------------------
% Saying you had some new resistance added
% could you differentiate them?
%-----------------------------------------
% come up with some new EoR
EoR_new = 3.3;

% Add the new EoR to the old data to shift it
newValues = zeros(1, (curveDataEnd-curveStart+1));
for i = (1:length(times))
    newValues(i) = flow(i+curveStart) * exp(times(i)*EoR_new);
end

% mess the new data up - pad the start and add some noise
padLength = 12;
padding = ones(1, padLength)*newValues(1);
messyValues = awgn([padding, newValues], 18);

% filter it a bit
sigma = 4;
size = 15;
x = linspace(0, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
filteredValues = filter(gaussFilter, 1, messyValues);

% remove padding from start
filteredValues = filteredValues(padLength:end);

figure(1)
hold on
plot(filteredValues, 'r', 'linewidth', 2)
legend("original", "LSQ fit", "extra R");
hold off
