% One thing to note here is that the file will have a colum Marker
% It will always have 3.something

% Video Group: 
% Love: 3.01, 3.02, 3.03
% Sad: 3.04, 3.05, 3.06
% Fear: 3.07, 3.08, 3.09
% Frustration: 3.10,3.11,3.12
% Calm: 3.13,3.14,3.15


    %% Variables setup
    

%% Load the raw data (From a folder) and iterate through each participant
%  For now we will only load one data file and play with it

    current_participant_id = 1;
    %% Calculate the features and return a matrix
    X = calculate_features(raw_data);
    sample_id = repmat(current_participant_id,[length(X) 1])
    
    Y = [Y;sample_id];
    
    %% Run the K-means
    idx = kmeans(X,k)
    % idx will contain the cluster belongings to a given sample!
    
    %% Assign a cluster to a participant depending on its ID

