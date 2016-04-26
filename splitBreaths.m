% Author: Sarah Howe
% Written in octave 4.01
% May not be compatible with matlab

function indices = splitBreaths(data)
    START = 1;
    END = 2;
    MIN_PEAK_VALUE = 0;
    if max(data(1:1000)) > 0.5 %flow data
        MIN_PEAK_VALUE = 0.15;
    else % pressure data
        MIN_PEAK_VALUE = 0.015;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % separate breaths with start of
    % inhalation at start, and end of
    % expiration at the end. Return
    % list of start and end indices for
    % each breath
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % A small half breath is about 1 s,
    % sampling frequency is 125 Hz, so guess
    % about 125 data points per half breath.
    % Take a percentage of that as test range
    % for start/end of breath
    data_points_per_breath = 125;
    percentage_of_breath_to_agree = 0.03;
    test_range = ceil(data_points_per_breath * percentage_of_breath_to_agree);

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
    for current_index = test_range+1:length(data)
        % Decide if we're looking for the
        % start or the end of a breath
        if(start_found)
            part_of_breath = END;
        else
            part_of_breath = START;
        end
        crossing_found = 0;

        % Look for crossing points or high data values
        if(part_of_breath == START);
            % look for a positive-gradient zero-crossing
            crossing_found = (data(current_index - test_range) <= 0 && data(current_index) > 0);
        else
            % look for a negative-gradient zero-crossing
            crossing_found = (data(current_index - test_range) >= 0 && data(current_index) < 0);
        end

        if(crossing_found);
            if(breath_count!= 1);
                % Check it has been a while since the last crossing
                if(current_index - breaths(part_of_breath, breath_count-1) > test_range);
                    breaths(part_of_breath, breath_count) = current_index;
                    start_found = ! start_found;
                    if(part_of_breath == END);
                        breath_count = breath_count + 1;
                    end
                end
            else
                % If it is the first breath, just save it
                breaths(part_of_breath, breath_count) = current_index;
                start_found = ! start_found;
                if(part_of_breath == END);
                    breath_count = breath_count + 1;
                end
            end
        end
    end

    % Remove any pairs which don't have a peak inbetween
    % or are too close together
    tempIndices = zeros(2, breath_count);
    position = 1;
    for pair = 1:breath_count
        % Get the start and end of each pair
        start = breaths(START, pair);
        stop = breaths(END, pair);
        % There has to be an end for the pair
        if((stop != 0) && (stop - start > test_range))
            for i = start:stop
                % If there is a peak point
                if abs(data(i)) > MIN_PEAK_VALUE
                    % Save the pair
                    tempIndices(:, position) = [start; stop];
                    position = position + 1;
                    break;
                end
            end
        end
    end
    % Cut off trailing zeros
    indices = tempIndices(:, 1:position-1);
end
