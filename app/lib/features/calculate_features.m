function [X] = calculate_features(raw_data)
%CALUCLATE_FEATURES Summary of this function goes here
%   Detailed explanation goes here

    %% Variables Declaration
    mean_window_size = 4;
    mean_window_overlap = 2;

    %% Blood Volume Pulse Features
    bvp_mean = running_mean("bvp",raw_data,mean_window_size,mean_window_overlap);
    %% Skin Conductance Features
    sc_mean = running_mean("sc",raw_data,mean_window_size,mean_window_overlap);
    %% Temperature Features
    temp_mean = running_mean("temp",raw_data,mean_window_size,mean_window_overlap);
    %% Heart Rate Features
    hr_mean = running_mean("hr",raw_data,mean_window_size,mean_window_overlap);
    
    X = [bvp_mean;sc_mean;temp_mean;hr_mean]';
end

