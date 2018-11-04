function [running_features] = iterate_over_windows(type,data_matrix,window_size,overlap,analysis_techniques)
%ITERATE_OVER_WINDOWS Summary of this function goes here
%   Input:
%       type: bvp,hr,sc,temp
%       data_matrix: the full dataframe
%       windows_size: size on which we want to calculate the analysis
%       techniques
%       overlap: overlap between each window
%       analysis_techniques: Structure containing a boolean for each
%       analysis techniques we want to implement (1 = do it, 0 = don't)
%   Output:
%       running_feature: structure containing vectors of running features

    %% Setting up the return structure
    running_features = struct();
    
    %% Unpacking the analysis techniques
    is_mean = analysis_techniques.is_mean;
    is_std = analysis_techniques.is_std;
    
    %% Setting up the variables
    index_type = get_index_type(type);    
    data = data_matrix(:,index_type);
    time_length = get_length_data(data_matrix);
    data_length = length(data);
    number_window = floor(time_length/window_size);
    points_per_window = floor(data_length/number_window);
    points_per_overlap = floor(data_length/(number_window*overlap));
    number_window = number_window*overlap;
    
    %% Iterate over the windows
    for i = 1:number_window-1
        % Get the right index for the window
        start_index = (i-1)*points_per_overlap + 1;
        end_index = (i-1)*points_per_overlap + points_per_window;
        
        % Analysis techniques:
        % Mean
        if(is_mean) 
            running_avg(i) = mean(data(start_index:end_index));
        end
        
        % Std
        if(is_std)
           running_std(i) = std(data(start_index:end_index));s 
        end
        
    end
    %% TODO check what to do with the remainder of the data (as we floor down)
    
    
    %% Load the returning stucture
    running_features.mean = running_avg;
    running_features.std = running_std;
    
end

