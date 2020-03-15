clear

datafile = load('SpirometryData.mat');
data = datafile.data;

% Has data for three people
people_to_keep = [1,2,3];

% Has sub sets of data
sets_to_keep = {'Loops','Banding','Normal','Inflated'};

% Has data in sets
data_to_keep = {'Pressure', 'Flow'};

for i = 1:length(sets_to_keep)
    for j = 1:length(people_to_keep)
        filename = [sets_to_keep{i} '_' int2str(people_to_keep(j)) '.csv']

        for dataindex = 1:length(data_to_keep)
            set = sets_to_keep{i};
            guy= people_to_keep(j);
            datatype = data_to_keep{dataindex};

            writedata = data(guy).(set);
            writedata = writedata.(datatype);
            dataindex;
            data_to_write(dataindex,:) = writedata;
        end

        fid = fopen(filename, 'wt');

        [cols, rows] = size(data_to_write);

        fprintf(fid, '%s,', data_to_keep{1,1:end-1})
        fprintf(fid, '%s\n', data_to_keep{end})

        for z = 1:rows
            fprintf(fid, '%f,', data_to_write'(z,1:end-1))
            fprintf(fid, '%f\n', data_to_write'(z,end))
        end
        fclose(fid);
        data_to_write = [];
        printf('%s\n', 'file created')
    end
end
