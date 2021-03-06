function [participants_clusters,percentages_aggrements] = run_clustering(X,Y,type,distance,start_cluster,end_cluster)
%RUN_CLUSTERING Will iterate from start to end and run the clustering
%analysis on all the possible cluster
%   Input:
%       type: type of clustering analysis to run
%       X: contain the features (sample x number of features)
%       Y: contain the id of which samples belong to which participant
%       start/end_cluster: first and last cluster to try out
%   Output:
%       participants_cluster: matrix of all the clustering id
%       percentages_aggrements: matrix of all the percentage aggrements

    %% Setting up variables
    participants_clusters = [];
    percentages_aggrements = [];
    
    %% Select which clustering technique to run
    is_kmeans = strcmp(type,"kmeans");
    
    %% Repeat for start_cluster, start_cluster+1...until end_cluster
    for k = start_cluster:end_cluster

        if(is_kmeans)
            disp("Running k means on " + num2str(k) + " clusters"); 
            idx = kmeans(X,k);
        end
        
        % Assign a cluster to a participant depending on its ID
        [participants_cluster,percentages_aggrement] = assign_cluster(idx,Y);
        
        % Store data
        participants_clusters = [participants_clusters,participants_cluster];
        percentages_aggrements = [percentages_aggrements,percentages_aggrement];
    end
end

