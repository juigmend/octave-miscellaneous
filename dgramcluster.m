function c = dgramcluster(Z,p,perm)
% c = dgramcluster(Z,p,perm)
%
% dgramcluster groups observations into clusters.
%
% c is a m x 2 matrix with observations indexes in the first column and
% extracted clusters in the second column.
%
% Z is a hierarchical tree (e.g., made with the linkage function).
% p is a dendrogram (e.g., made with the dendrogram function).
% perm is the permutation of observations (e.g., made with the dendrogram function).
%
% Juan Ignacio Mendoza Garay - 2015

obs_amt = size(Z,1) + 1; % amount of observations
d_clustering = squareform(pdist(Z(:,3))); % distance matrix of linkage distances
first_diag = (diag(d_clustering,1))'; % 1st. diagonal's maxima are cluster separations
[~, tree_thresh] = max(first_diag); % 1st. diagonal's maximum is threshold for biggest clusters
det_clusts = obs_amt - tree_thresh; % amount of more relevant (biggest) clusters
leaves_n_clusts = p(1:(end-det_clusts+2)-1,1:2); % extract x position of leaves and clusters
bounds = leaves_n_clusts(diff(leaves_n_clusts(:,2))==0); % cluster boundaries are at two consecutive leaves
c(:,1) = perm;
c(:,2) = det_clusts;
counter = 1;
for i = 1:det_clusts-1 % assign labels to clusters
    c(counter:bounds(i),2) = i;
    counter = bounds(i) + 1;
end



