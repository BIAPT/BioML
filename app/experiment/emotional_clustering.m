% One thing to note here is that the file will have a colum Marker
% It will always have 3.something

% Video Group: 
% Love: 3.01, 3.02, 3.03
% Sad: 3.04, 3.05, 3.06
% Fear: 3.07, 3.08, 3.09
% Frustration: 3.10,3.11,3.12
% Calm: 3.13,3.14,3.15


    %% Variables setup;
    extension = "*.csv";
    X = [];
    Y = [];

    %% Load the raw data (From a folder) and iterate through each participant
    %  For now we will only load one data file and play with it
    data_folder = uigetdir;
    data_files = dir(fullfile(data_folder,extension)); %this is a structure
    
    % Iterating over the participant and calculating the features
    for participant_id = 1:length(data_files)
        %% Calculate the features and return a matrix
        file_name = strcat(data_files(participant_id).folder,filesep,data_files(participant_id).name);
        raw_data = load_data(file_name);
        features_matrix = calculate_features(raw_data);
        sample_id = repmat(participant_id,[length(features_matrix) 1]);
        
        X = [X;features_matrix];
        Y = [Y;sample_id];
    
    end
    
    %% Run the K-means
    [participants_clusters,percentages_aggrements] = run_clustering("kmeans",X,Y,4,34);
   