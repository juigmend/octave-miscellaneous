%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                    Supervised Classification of Time Series                  %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             Doctoral Student %
%                                   Music Department - University of Jyväskylä %
%                                                                   July, 2015 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;
% ------------------------------------------------------------------------------
% INPUT DATA AND PARAMETERS:

load 'labelled_series_train' % <-- load labelled training timeseries
load 'labelled_series_test'  % <-- load labelled testing timeseries
winperc = 10; % <----------------- DTW window is a percentage of timeseries length

% Notes on the input: 
% - All timeseries should have the same length in amount of samples (rows).
% - First column is labels.
% - Classes should be distinct enough and timeseries within each class should be
%   similar enough. How much is enough? That's what we want to find.

% Notes on the method:
% - Dynamic Time Warping is used as distance measure.
% - One Nearest Neighbour is used as classifier.
% - Other distance and classification methods could be used although these 
%   seem to work very well with this kind of data.

% Dependencies:
% - Octave (maybe also runs in Matlab).
% - An external Dynamic Time Warping function.

% ------------------------------------------------------------------------------

train_set = series_train;
test_set = series_test;
train_labels = train_set(:,1); % extract training class labels
train_set(:,1) = []; % remove class labels from training set
test_labels = test_set(:,1); % extract testing class labels
test_set(:,1) = []; % remove class labels from test set
w = round((size(train_set,2))*winperc/100); % size of DTW window
amt_test_labels = length(test_labels);
amt_train_labels = length(train_labels);
report = zeros(amt_test_labels,2);
report(:,1) = test_labels;

wbar = waitbar(0,'patience is a virtue');
fractot = 1/amt_test_labels; 
laps = 1;
for i = 1 : amt_test_labels % compare testing series with training series
   % --- 1NN classifier ---
   control = inf;
   for i_1 = 1:amt_train_labels
     d = dtw(train_set(i_1,:)', test_set(i,:)', w); % distance function
     if d < control
       detected_class = train_labels(i_1);
       control = d;
     end
   end;
   % --- 
   report(i,2) = detected_class;
   waitbar(fractot*laps);
   laps = laps+1;
end;
close (wbar)

% generate report:
report
fprintf('[original class, detected class] \n')
fprintf('Accuracy = %i%% \n',...
 (length(find((report(:,1)-report(:,2))==0)))*100/amt_test_labels)



