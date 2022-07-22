%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%            Unsupervised Non-Parametric Clustering of Time Series             %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             Doctoral Student %
%                                   Music Department - University of Jyväskylä %
%                                                                   July, 2015 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all
% ------------------------------------------------------------------------------
% INPUT DATA AND PARAMETERS:
 
load 'labelled_series' % <------ labelled randomly sorted n-classes timeseries
winperc = 10; % <--------------- DTW window is a percentage of timeseries length

% Notes on the input: 
% - All timeseries should have the same length in amount of samples (rows).
% - First column is labels.
% - Classes should be distinct enough and timeseries within each class should be
%   similar enough. How much is enough? That's what we want to find.

% Notes on the method:
% - Dynamic Time Warping is used as distance measure.
% - One Nearest Neighbour is used as hierarchical clustering.
% - A heuristic of maximum distance between nodes is used to identify clusters.
% - Other distance and clustering methods could be used although these seem to
%   work very well with this kind of data.

% Dependencies:
% - Octave (maybe also runs in Matlab).
% - Octave statistics package (maybe something similar is available for Matlab).
% - An external Dynamic Time Warping function.

% ------------------------------------------------------------------------------

class_labels(:,2) = series(:,1); % extract class labels
class_labels(:,1) = [1:size(class_labels,1)]; % index class labels
amt_classes = size(unique(class_labels(:,2)),1);
series(:,1) = []; % remove class labels from data matrix
series_length = size(series,2);
amt_ser = size(series,1);
w = round(series*winperc/100); % size of DTW window

% make DTW distance matrix of all signals:
wbar = waitbar(0,'patience is a virtue');
fractot = 2/(amt_ser^2 - amt_ser); 
laps = 1;
d = zeros(amt_ser);
d_half = zeros(amt_ser);

for row = 2:amt_ser
  for col = 1:(row-1)
    d(row,col) = dtw(series(row,:)',series(col,:)',w); % distance function
    d(col,row) = d(row,col);
    d_half(row,col) = d(row,col);
    waitbar(fractot*laps);
    laps = laps+1;
  end
end
close(wbar)

imagesc(d) % plot heat map of distance matrix
d_half_vec = reshape(d_half,1,amt_ser^2);
d_vec = nonzeros(d_half_vec);
Z = linkage(d_vec,'single'); % hierarchical clustering
figure
[p, ~, perm] = dendrogram(Z); % display dendrogram
hold on

% ------------------------------------------------------------------------------
% THIS IS MY HUMBLE CONTRIBUTION:
% A very simple method to determine the major clusters based in a heuristic 
% of maximum distance between clusters is used to identify clusters., which should be the 
% classes of the data. This will work only if the different classes are distinct
% enough and the series within each class are similar enough.
% You can choose if

% a) you want to use the function. Comment/Uncomment the following line.
%c = dgramcluster(Z,p,perm);

% b) you want to use the raw code (same as in the function)
% so you can inspect it. Comment/uncomment where indicated.
% Some following code is redundant but we need it here to make the beautiful
% plot over the dendrogram.

% extract amount of main clusters:
d_clustering = squareform(pdist(Z(:,3))); % distance matrix of linkage distances
first_diag = (diag(d_clustering,1))'; % 1st. diagonal's maxima are cluster separations
[~, tree_thresh] = max(first_diag); % 1st. diagonal's maximum is threshold for biggest clusters
det_clusts = amt_ser - tree_thresh; % amount of more main clusters 

% find elements of clusters:
leaves = p(1:amt_ser,1:2); 
clusters_nodes = p(amt_ser+1:(end-det_clusts+2)-1,1:2); 
classes_nodes = p((end-det_clusts+2):end-1,1:2);  
universe_node = p(end,1:2); % 
leaves_n_clusts = sortrows(cat(1,leaves,clusters_nodes),1); % boundaries are at duplicated zeroes

% Comment/uncomment the following:
bounds = leaves_n_clusts(diff(leaves_n_clusts(:,2))==0); % cluster boundaries
c(:,1) = perm;
c(:,2) = det_clusts;
counter = 1;
for i = 1:det_clusts-1 % assign labels to clusters
    c(counter:bounds(i),2) = i;
    counter = bounds(i) + 1;
end

% ------------------------------------------------------------------------------

% plot nodes superimposed to dendrogram:
plot(universe_node(1),universe_node(2),'k.','markersize',20);
plot(classes_nodes(:,1),classes_nodes(:,2),'g.','markersize',20);
plot(clusters_nodes(:,1),clusters_nodes(:,2),'r.','markersize',20);
plot(leaves(:,1),leaves(:,2),'b.','markersize',20);

% make parabola for animated plot of clusters:
x_big_circ = [-1:0.1:1.1]; % x
y_big_circ = sqrt(1-x_big_circ.^2) * clusters_nodes(1,2)/2; % y
counter = 1;
binary = 1;
for i = 1:size(y_big_circ,2) % shrink y to half
    if binary == 1
      y_small_circ(counter) = y_big_circ(i);
      counter = counter + 1;
      binary = 0;
    else
      binary = 1;
    end
end
x_small_circ = [0:0.1:1]; % shrink x

% animated plot of clusters:
% By doing this plot I discovered that there are leaves that are not joined,
% which appear as two consecutive zero-distances in the dendogram when sorted
% by the x values (first column). The first version of the code showed a 
% horizontal line between the cluster boundaries. Then I made the thing to
% skip the cluster boundaries and to make the nice parabolas :-)
for i = 1:size(leaves_n_clusts,1)-1
    if leaves_n_clusts(i,2)== 0 
      if leaves_n_clusts(i+1,2) ~= 0
        plot(x_small_circ + leaves_n_clusts(i,1),...
        y_small_circ,'k--','linewidth',6)
        pause(0.5)
      end
    end
end
hold off

% generate report:
report = zeros(amt_ser,4);
report = cat(2,c,sortrows(class_labels,2))
fprintf('[index, detected class, index, original class] \n')
fprintf('Detected and original class labels might be different as the algorithm \n')
fprintf('was not shown the original class labels. What matters is the grouping ;-) \n')
fprintf('Detected amount of classes: %i \n',det_clusts)
fprintf('True amount of classes: %i \n',amt_classes)





