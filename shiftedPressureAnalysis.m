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
Fs = 125;

% start loops
flow = loops.Flow;
pressure = loops.Pressure;
%flow = normal.Flow;
%pressure = normal.Pressure;

% unit conversions
L_to_ml = 1000;
s_to_min = 1/60;
kPa_to_cmH20 = 10.1972;

%Pressure: kPa -> cmH2O
pressure = pressure*kPa_to_cmH20;

% time for plotting
time = (1:size(flow))*(1/Fs);

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

%{
% plot pressure and reduced flow with split points
figure(1)
hold on
plot(time, -pressure, '.-b')
plot(time, flow/10, '.-k')
for value = 1:length(pressure_splits)
    plot(time(pressure_splits(START, value)), -pressure(pressure_splits(START, value)), '.g')
    plot(time(pressure_splits(END, value)), -pressure(pressure_splits(END, value)), '.r')
end
for value = 1:length(flow_splits)
    plot(time(flow_splits(START, value)), flow(flow_splits(START, value))/10, '.c')
    plot(time(flow_splits(END, value)), flow(flow_splits(END, value))/10, '.m')
end
legend('-Pressure', 'Flow x 0.1')
grid minor
hold off
%}

% Start of insp and exp stored in matrix
INHALE = 1;
EXHALE = 2;

% Get first pair from pressure and flow
pressure_pair = 1;
flow_pair = 1;

% start and end indices of inspiration
%flow_start = flow_splits(INHALE, flow_pair);
%flow_end = flow_splits(EXHALE, flow_pair);

% Look at all of the breath pairs in the data
storage_index = 1;
breath_state = INHALE;
pairs = zeros(2, length(pressure_splits));
while(pressure_pair < length(pressure_splits) && flow_pair < length(flow_splits))
    % start and end indices of inspiration
    pressure_start = pressure_splits(INHALE, pressure_pair);
    pressure_end = pressure_splits(EXHALE, pressure_pair);

    % Work out length of whole breath assuming
    % inh ans exh are the same length and shape
    % for a breath
    breath_length = 2*(pressure_end - pressure_start);

    % Check that the flow and pressure begin within
    % a half breath of each other (with compliance it
    % can't be further than a half breath out)
    % Function returns corrected indices if breath
    % missing in either set.
    [pressure_pair_corrected, flow_pair_corrected] = checkIndicesAreForSameBreath(pressure_pair, ...
                            pressure_splits, flow_pair, flow_splits, breath_length, breath_state);

    % The function will return non-indexed values
    % if no pair was found.
    if(pressure_pair_corrected > length(pressure_splits) || flow_pair_corrected > length(flow_splits));
        disp('No matching index pair found for breath');
        break;
    end

    % record correct pairs
    pairs(:, storage_index) = [pressure_pair_corrected; flow_pair_corrected];

    % Go to the next half breath
    if(breath_state == EXHALE);
        breath_state = INHALE;
        pressure_pair = pressure_pair + 1;
        flow_pair = flow_pair + 1;
    else
        breath_state = EXHALE;
    end

    storage_index = storage_index + 1;
end

% Drop the trailing zero pairs
pairs = pairs(:,1:storage_index-1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Find phase shift from fourier transform
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
breath_state = INHALE;
for breath = 1:length(pairs)
    %range being looked at is a half breath
    if(breath_state == INHALE)
        % start and end indices of inspiration
        pressure_start = pressure_splits(INHALE, breath);
        pressure_end = pressure_splits(EXHALE, breath);
        range = (pressure_start:pressure_end;
    else
        if(breath < length(pairs))
            % start and end indices of inspiration
            pressure_start = pressure_splits(EXHALE, breath);
            pressure_end = pressure_splits(INHALE, breath+1);
            range = (pressure_splits(2,breath):pressure_splits(1,breath+1));
        else
            
        end
    end

    % Will add inverted half breath on the end to make a full waveform
    length_range = length(range)*2;

    % Frequency for plotting
    freq = Fs*(0:length_range/2)/length_range;

    % find the data point difference between flow and pressure
    dif = (flow_splits(breath_state, pairs(2,breath))-pressure_splits(breath_state, pairs(1,breath)))

    % Shift breath to frequency domain
    % Get half the breath since mirrored
    padding = zeros(dif, 1);
    flow_range = range + ones(1,length(range))*dif;
    flow_curve = [padding; flow(flow_range); -flow(flow_range)];
    pressure_curve = [-pressure(range); pressure(range); padding];

    % fourier transfrom
    fftF = fft(flow_curve);
    fftP = fft(pressure_curve);

    % get the normalised absolute values
    freqF = abs(fftF/length_range);
    freqP = abs(fftP/length_range);

    % take the first half, since repeated
    freqF = freqF(1:length_range/2+1);
    freqP = freqP(1:length_range/2+1);

    % double magnitude to make up for lost points
    freqF(2:end-1) = (2*freqF(2:end-1));
    freqP(2:end-1) = (2*freqP(2:end-1));


    %{
    figure(2)
    subplot(211)
    plot(freq, freqF)
    ylabel('Flow (L/s)')
    xlabel('Frequency Hz')
    grid minor
    subplot(212)
    plot(freq, freqP)
    ylabel('Pressure (cmH20)')
    xlabel('Frequency Hz')
    grid minor
    %}

    % Find the fundamental frequency
    [magnitude_flow, max_index_flow] = max(freqF);
    [magnitude_pressure, max_index_pressure] = max(freqP);

    % Find the phases
    phase_flow = angle(fftF(max_index_flow));
    phase_pressure = angle(fftP(max_index_pressure));
    phase_shift = phase_pressure - phase_flow
    phase_shift_degrees = phase_shift*180/pi

    % shift the flow so it aligns with pressure
    shifted_range = range(dif+1:end).-ones(1,length(range)-dif)*dif;

    % plot the weird flipped data with the phase estimates
    np = zeros(1, length_range);
    np(max_index_pressure) = fftP(max_index_pressure);
    npt = ifft(np);
    nf = zeros(1, length_range);
    nf(max_index_flow) = fftF(max_index_flow);
    nft = ifft(nf);
    figure(3)
    hold on
    plot([-pressure(range);pressure(range)], 'b');
    plot(real(npt), '--b')
    plot([flow(range)/1;-flow(range)/1], 'k');
    plot(real(nft)/1, '--k')
    grid minor
    legend('pressure data', 'pressure phase', 'flow data', 'flow phase')
    hold off

    % Real part of complex reactance
    % remove data point difference from start of flow
    % remove data point difference from end of pressure
    % complex_resistance = V/I


    % find the resistance using the aligned data
    resistance = flow(shifted_range)\-pressure(range(dif+1:end))

    % resistance and reactance separated from complex value
    %resistance = complex_resistance*cos(phase_shift)
    reactance = resistance*tan(phase_shift)

    % Get the compliance from reactance
    signal_freq = Fs/length_range;
    elastance = -reactance*(2*pi*signal_freq)

    e(breath) = elastance;
    r(breath) = resistance;

    if(breath_state == INHALE)
        breath_state = EXHALE
    else
        breath_state = INHALE
    end

    figure(5);
    hold on
    plot(pressure_curve, '.-')
    plot(flow_curve, '.-r')
    title(breath)
    legend('pressure','flow')
    grid minor
    hold off

    pause()

    close(5)
    close(3)
end

figure(6)
hold on
plot(1:breath, -e, 'ob')
plot(1:breath, r, 'xr')
xlabel('breath')
ylabel('magnitude')
legend('elast', 'resist')
grid minor
hold off

%{
r = 1:1000;
p = flow(r)*resistance;
P = fft(p);
shift_as_ratio = phase_shift/(2*pi);
f=linspace(0, Fs, r);
P = P.*exp(-j*2*pi*f*shift_as_ratio);
calcP = ifft(P);

figure()
hold on
plot(pressure(r), 'b')
plot(real(calcP), 'r')
plot(flow(r)*resistance, 'k', 'linewidth', 3)
legend('old pressure', 'calculated pressure', 'flow')
grid minor
hold off
%}
