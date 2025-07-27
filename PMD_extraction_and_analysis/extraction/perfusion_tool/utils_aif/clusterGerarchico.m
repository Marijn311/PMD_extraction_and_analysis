
%% ------------------------------------------------------------------------
function [clusterAssignments, centroids] = clusterGerarchico(data,nClusters)
    % Applies hierarchical clustering to divide data into nClusters
    %
    % Inputs:
    %   data - Matrix where each row is a time series to cluster
    %   nClusters - Number of clusters to create
    %
    % Outputs:
    %   clusterAssignments - Vector indicating cluster assignment for each row
    %   centroids - Matrix containing mean time series for each cluster
    
    % Calculate pairwise distances and create hierarchical cluster tree
    distance = pdist(data);
    tree = linkage(distance,'ward');
    
    % Cut tree to get desired number of clusters
    clusterAssignments = cluster(tree,'maxclust',nClusters);
    
    % Calculate centroids (mean time series) for each cluster
    nT = size(data,2);
    centroids = zeros(nClusters,nT);
    for k=1:nClusters
        ind = find(clusterAssignments==k);
        
        clusterData = zeros(length(ind),nT);
        for t=1:nT
            clusterData(:,t) = data(ind,t);
        end
        
        centroids(k,:) = mean(clusterData,1);
    end
    end % clusterGerarchico
    
    