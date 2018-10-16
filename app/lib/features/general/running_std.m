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
end

