function print_clusters_ratio(clusters)
%GET_CLUSTERS_RATIO Summary of this function goes here
%   Detailed explanation goes here

    %% Variables Declaration
    
    unique_clusters = unique(clusters);
    proportion_clusters = [];
    for i = 1:length(unique_clusters)
       number_current_cluster = length(clusters(clusters==unique_clusters(i)));
       proportion_clusters = [proportion_clusters;number_current_cluster];
    end
    
    figure;
    title("Participants Cluster Ratio");
    pie(proportion_clusters,cellstr(num2str(unique_clusters)));
    
end

