% Author: Sarah Howe
% Written in octave 4.01
% May not be compatible with matlab

function [pressure_pair, flow_pair] = checkIndicesAreForSameBreath(pressure_pair, ...
                                                                pressure_splits, ...
                                                                flow_pair, ...
                                                                flow_splits, ...
                                                                breath_length)
    noMatch = 1;
    while(noMatch)
        % Find indices of starting points
        pressure_start = pressure_splits(1, pressure_pair);
        flow_start = flow_splits(1, flow_pair);
        % Check if they are close enough
        if(abs(pressure_start - flow_start) < breath_length/2)
            noMatch = 0;
        else
            % Go to next pressure index if flow leads
            if(pressure_start < flow_start)
                pressure_pair = pressure_pair + 1;
                if(pressure_pair > length(pressure_splits))
                    noMatch = 0;
                end
            else
                % Go to next flow index
                flow_pair = flow_pair + 1;
                if(flow_pair > length(flow_splits))
                    noMatch = 0;
                end
            end
        end
    end
end
