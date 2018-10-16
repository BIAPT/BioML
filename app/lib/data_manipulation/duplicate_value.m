function [replicate_vector] = duplicate_value(value_to_duplicate,number_duplication)
%DUPLICATE_VALUE Wrapper function to duplicate a value
%   Input:
%       value_to_duplicate: this is the value we wish to duplicate
%       number_duplication: this is the number of time we want to duplicate
%       it
%   Output:
%       replicate_vector: the replicated value

   replicate_vector = repmat(value_to_duplicate,1,number_duplication);
end

