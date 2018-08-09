function [data_struct] = import_csv(file_path)
%IMPORT_CSV Import a csv files containing biosignal data into a table
%   Detailed explanation goes here

% TODO create a data structure that will hold the data for each condition
% data_struct. Here we want to keep as much information as possible and not
% process too much.
%   0 (start)
%   1 (collect_qc_data)
%   2 (run_qc)
%   3 (video_playing)
%   4 (video_testing)
%   5 (instructions)
%   6 (replay)


data_table = readtable(file_path);

end

