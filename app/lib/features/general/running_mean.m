function [running_avg] = running_mean(type,data_matrix,length_mean,overlap)
%RUNNING_MEAN Summary of this function goes here
%   Detailed explanation goes here
    
    index_type = 0;
    if(strcmp(type,"bvp"))
        index_type = 2;
    elseif(strcmp(type,"sc"))
        index_type = 3;
    elseif(strcmp(type,"temp"))
        index_type = 4;
    elseif(strcmp(type,"hr"))
        index_type = 5;
    end
    time_length = get_length(data_matrix);
    
end

