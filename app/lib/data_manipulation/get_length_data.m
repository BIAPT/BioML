function [time_length] = get_length_data(counter)
%GET_LENGTH_DATA Summary of this function goes here
%   Detailed explanation goes here
%% Need to return the length in seconds
    %May need to convert
    time_length = counter(length(counter),1) - counter(1,1);

end

