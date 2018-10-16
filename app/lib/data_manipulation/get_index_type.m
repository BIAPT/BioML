function [index_type] = get_index_type(type)
%GET_INDEX_TYPE Helper function to get the index in the column of a given
%type
%   Input:
%       type: bvp,sc,temp or hr
%   Output:
%       index_type: columns

    if(strcmp(type,"bvp"))
        index_type = 2;
    elseif(strcmp(type,"sc"))
        index_type = 3;
    elseif(strcmp(type,"temp"))
        index_type = 4;
    elseif(strcmp(type,"hr"))
        index_type = 5;
    end
end

