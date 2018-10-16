function [running_std] = running_std(type,data_matrix,length_std,overlap)
%RUNNING_STD calculatie the running std for a given modality
%   Input:
%       type: bvp,hr,sc,temp
%       data_matrix: the full dataframe
%       length_std: size on which we want to calculate the means
%       overlap: overlap between each window
%   Output:
%       running_std: the vector containing the running std

    index_type = get_index_type(type);
    data = data_matrix(:,index_type);
    time_length = get_length_data(data_matrix);
    data_length = length(data);
    number_window = floor(time_length/length_std);
    points_per_window = floor(data_length/number_window);
    points_per_overlap = floor(data_length/(number_window*overlap));
    number_window = number_window*overlap;
    for i = 1:number_window-1
        start_index = (i-1)*points_per_overlap + 1;
        end_index = (i-1)*points_per_overlap + points_per_window;
        running_std(i) = std(data(start_index:end_index));
    end
end

