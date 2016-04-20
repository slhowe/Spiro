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

% start loops
flow = loops.Flow;
pressure = loops.Pressure;

% sampling frequency 125 Hz
Hz = 125;

% time for plotting
time = (1:size(flow))*(1/Hz);

% Calc volume.
% This drifts over many breaths
% best to go breath to breath
vol = cumtrapz(flow*(1/Hz));

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

% Looks like an RC curve at end of forced exp?
% Cause? Lung relaxing after muscles release???

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

% set stopping point of fit at 90% flow drop
drop = 0.9*flow(curveStart)-flow(curveDataEnd);

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
curveStop = CurveStart+index;

% set up matrices
measurements = log(-flow(curveStart:curveStop));
one = ones(1, curveStop-curveStart+1);
times = -(time(curveStart:curveStop)-time(curveStart));

% OMG least squares!!!
results = [one', times']\measurements;

% extract info
startPoint = exp(results(1));
EoR = results(2);

% remake curve from info
newValues = startPoint*exp(times*EoR);

% compare new and old
figure(1)
hold on
plot(-newValues, 'm')
legend("original", "LSQ fit");
hold off


%E = zeros(1,1200);
%R = zeros(1,1200);
%pi = p_atm;
%for i = 4:1200
%    vi = vol(i) - offset;
%    if vi == 0
%        E(i) = 0;
%    else
%        E(i) = pi/vi;
%    endif
%    params = [vol(i-3:i)-offset, flow(i-3:i)]\(pi*ones(4, 1));
%    E(i) = params(1);
%    R(i) = params(2);
%endfor

%figure
%hold on
%subplot(4, 1, 1)
%hold on
%plot(time(1:1200), flow(1:1200))
%plot(time(740), flow(740), 'rx')
%hold off
%grid minor
%legend("flow")
%ylabel("l/s")
%xlabel("s")
%subplot(4, 1, 2)
%hold on
%plot(time(1:1200), vol(1:1200)-offset+3)
%plot(time(740), vol(740)-offset, 'rx')
%hold off
%legend("volume")
%grid minor
%ylabel("L")
%xlabel("s")
%subplot(4, 1, 3)
%hold on
%plot(time(1:1200), E)
%plot(time(740), E(740), 'rx')
%hold off
%grid minor
%ylabel("E")
%subplot(4, 1, 4)
%hold on
%plot(time(1:1200), R)
%plot(time(740), R(740), 'rx')
%hold off
%grid minor
%ylabel("R")

