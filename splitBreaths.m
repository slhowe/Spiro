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

function indices = splitBreaths(data)
    START = 1;
    END = 2;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % separate breaths with start of
    % inhalation at start and end of
    % expiration at the end. Return
    % list of start and end indices for
    % each breath
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Each half of a breath is about 1 s
    % sampling frequency is 125 Hz so
    % about 125 data points per half breath
    % take about 5% of that as test range
    % for start/end of breath
    sampling_freq = 125;
    percentage_of_breath_to_agree = 0.05;
    test_range = ceil(sampling_freq * percentage_of_breath_to_agree);

    % Can't have more breaths than 1/2
    % the number of the data points.
    % Have first row as start indices,
    % second row as end indices
    breaths = zeros(2, length(data)/2);

    % Will need to know the actual number
    % of breaths for later
    breath_count = 1;

    % Find breath start and end indices.
    start_found = 0;
    for i = test_range:length(data)
        % Decide if we're looking for the
        % start or the end of a breath
        if(start_found)
            part_of_breath = END;
        else
            part_of_breath = START; 
        end
        % look for zero crossing over range
        if(data(i - test_range) <= 0 && data(i) > 0)
            % Check it's not noise around the crossing
            if(i - breaths(part_of_breath, breath_count) > test_range
                breaths(part_of_breath, breath_count) = i;
                if (part_of_breath == END)
                    breath_count = breath_count + 1;
                end
                start_found = !start_found; % Toggle part of breath being stored
            end
        end
    end
end
