pkg load communications

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

%-----------------------------------------------
% Looks like an EoR curve at end of forced exp?
% Let's fit a line to it :)
%-----------------------------------------------

% EoR curve range
start = 1190;
stop = 2000;

figure(1)
hold on
plot(flow(start:stop), 'b', 'linewidth', 3)
xlabel("dataPoint")
ylabel("flow")
grid minor
hold off

% define start and end point of whole curve
curveStart = start;
curveDataEnd = stop;

% set stopping point of fit at certain %age flow drop
% (things get over constrained below 3% and fit gets bad)
drop = flow(curveDataEnd) + 0.1*(flow(curveStart)-flow(curveDataEnd));

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
one = ones(1, (curveStop-curveStart)+1);
times = -(time(curveStart:curveStop)-time(curveStart));

% OMG least squares!!!
results = [one', times']\measurements;

clc % had to add this again, so annoying =/

% extract info
startPoint = exp(results(1));
EoR = results(2)

% remake curve from info
times = -(time(curveStart:curveDataEnd)-time(curveStart));
newValues = -startPoint*exp(times*EoR); % flip flow to negative side

% compare new and old
figure(1)
hold on
plot(newValues, 'm',  'linewidth', 2)
hold off

%-----------------------------------------
% Saying you had some new resistance added
% could you differentiate them?
%-----------------------------------------
% come up with some new EoR
EoR_new = 20

% Add the new EoR to the old data to shift it
newValues = zeros(1, ((curveDataEnd-curveStart)+1));
for i = (1:length(times))
    newValues(i) = flow(i+curveStart) * exp(times(i)*EoR_new);
end

% mess the new data up - add some noise
messyValues = awgn(newValues, 28);

% remove any noise that crosses or lies on the x-axis
filteredValues = zeros(1, length(messyValues));
if messyValues(1) >= 0
    error("First data point is not usable - try different start point")
end
for i = 1:length(messyValues)
    if messyValues(i) >= 0
        filteredValues(i) = filteredValues(i-1);
    else
        filteredValues(i) = messyValues(i);
    end
end

figure(1)
hold on
plot(filteredValues, 'r', 'linewidth', 3)
hold off

%------------------------------------------------------
% Now we have some data with added resistance in series
% Try to separate it all again
%------------------------------------------------------
% have Q = Ae^(k*t)
% where A = Q(1)e^(EoR*t)

% when modelled flow goes below -1, the reduced flow explodes
% so reduced flow calculation is reduced to a smaller range
times = -(time(curveStart:curveStop)-time(curveStart));
filteredValues = filteredValues(1:(curveStop-curveStart)+1);

% model flow drop across spirometer's filter from known EoR
lowRFlow = -startPoint*exp(times*EoR);

% remove the known flow drop from flow
reducedFlow = filteredValues./lowRFlow;

figure(1)
hold on
plot(-reducedFlow, 'k', 'linewidth', 2)
plot(-exp(times*EoR_new), 'c', 'linewidth', 2)
hold off

% set stopping point of fit at certain %age flow drop
% (things get over constrained below 3% and fit gets bad)
curveLength = length(reducedFlow);
drop = reducedFlow(curveLength) + 0.05*(reducedFlow(1)-reducedFlow(curveLength));

% find the index of the stopping point for the new curve
index = 0;
stillLooking = 1;
for i = 1:curveLength
    if(stillLooking)
        if(reducedFlow(i) < drop)
            index = i;
            stillLooking = 0;
        end
    end
end
if index == 0
    error("Percentage flow drop specified not found in range")
end
curveStop = index;

% redefine time for percentage pressure drop
times = -(time(curveStart:(curveStart+curveStop))-time(curveStart));
reducedFlow = reducedFlow(1:curveStop+1);

% find EoR for added flow drop
EoR_est = times'\log(reducedFlow)'

figure(1)
hold on
plot(-exp(times*EoR_est), 'g', 'linewidth', 2)
legend("Flow data", ...
        "LSQ fit", ...
        "Extra R added in series", ...
        "Flow drop across added R", ...
        "Actual EoR added", ...
        "Calculated EoR added");
hold off

EoR_Percent_error = abs((1 - EoR_est/EoR_new)) * 100

%----------------------------
% Say added R = 1 cmH20
%----------------------------
R_addedInSeries = 1
E = EoR_est*R_addedInSeries
R_spirometerFilter = 1/(EoR/E)
