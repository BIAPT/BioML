function [data_matrix] = load_data(file_path)
%LOAD_DATA Summary of this function goes here
%   Detailed explanation goes here

    data_matrix = csvread(file_path,1,0);
end

