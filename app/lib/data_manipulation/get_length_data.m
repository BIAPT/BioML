function [length] = get_length_data(counter)
%GET_LENGTH_DATA Summary of this function goes here
%   Detailed explanation goes here
%% Need to return the length in seconds
    %May need to convert
    length = counter(1) - counter(length(counter));

end

