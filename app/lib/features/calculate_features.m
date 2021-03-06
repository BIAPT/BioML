function [X] = calculate_features(raw_data)
%CALCuLATE_FEATURES will calculate all the features in a loop
%   Input:
%       raw_data: the unprocessed raw data
%   Output:
%       X: Features matrix

    %% Variables Declaration
    window_size = 4;
    window_overlap = 2;

    %% Blood Volume Pulse Features
    bvp_mean = running_mean("bvp",raw_data,window_size,window_overlap);
    bvp_std = running_std("bvp",raw_data,window_size,window_overlap);
    %% Skin Conductance Features
    sc_mean = running_mean("sc",raw_data,window_size,window_overlap);
    sc_std = running_std("sc",raw_data,window_size,window_overlap);
    %% Temperature Features
    temp_mean = running_mean("temp",raw_data,window_size,window_overlap);
    temp_std = running_std("temp",raw_data,window_size,window_overlap);
    %% Heart Rate Features
    hr_mean = running_mean("hr",raw_data,window_size,window_overlap);
    hr_std = running_std("hr",raw_data,window_size,window_overlap);
    
    X = [bvp_mean;bvp_std;sc_mean;sc_std;temp_mean;temp_std;hr_mean;hr_std]';
end

