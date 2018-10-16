function [running_avg] = running_mean(type,data_matrix,length_mean,overlap)
%RUNNING_MEAN Summary of this function goes here
%   Detailed explanation goes here
    
    index_type = get_index_type(type);    
    
    data = data_matrix(:,index_type);
    time_length = get_length_data(data_matrix);
    data_length = length(data);
    number_window = floor(time_length/length_mean);
    points_per_window = floor(data_length/number_window);
    points_per_overlap = floor(data_length/(number_window*overlap));
    number_window = number_window*overlap;
    for i = 1:number_window-1
        start_index = (i-1)*points_per_overlap + 1;
        end_index = (i-1)*points_per_overlap + points_per_window;
        running_avg(i) = mean(data(start_index:end_index));
    end
    
    %% TODO check what to do with the remainder of the data (as we floor down)
    
end

