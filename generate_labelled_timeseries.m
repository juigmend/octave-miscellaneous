%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                             Generate and Plot                                %
%                   N-classes Labelled Eigenmodes Time Series                  %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             Doctoral Student %
%                                   Music Department - University of Jyväskylä %
%                                                                   July, 2015 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
% ------------------------------------------------------------------------------
% Generate artificial time series that mimick bodily movement to music with 
% clear beat.
% ------------------------------------------------------------------------------
% INPUT VARIABLES:

beats = 8;  % <---------- amount of beats
srate = 10; % <---------- sampling rate in Hz
periods = [2 4 8]; % <--- period of the classes (1/Eigenfrequencies)
amt_ser = [2 2 2]; % <--- amount of series for each class
reorder_rand = 1;  % <--- reorder time series randomly (1 = yes, 0 = no)
dataset_label = 1; % <--- label of the dataset
% 0 = series
% 1 = series_train
% 2 = series_test

% Note: labels will be the indexes of "periods", same as "amt_ser".
% ------------------------------------------------------------------------------

tot_amt_ser = sum(amt_ser);
amt_classes = size(amt_ser,2);

x = [0:(1/srate):beats]; % time grid 

close all
figure
set(gca, 'XTick',[0:1:16]); grid on; hold on
cm = colormap(hsv(64));
cm_step = round(size(cm,1)/amt_classes);

% make and plot signals:
series = zeros(tot_amt_ser,size(x,2)+1);
counter = 1;
for i_classes = 1:amt_classes % classes
  colour = cm(1 + (i_classes-1)*cm_step,:);
    for i_series = 1:amt_ser(i_classes) % series per class

      series(counter,1) = i_classes; % first row has the index
    
      series(counter,2:end) = sin(x*pi*2/(periods(i_classes)) + randn(1)/8)...
      +(sin(x*pi/4*randn(1)*2)/6+randn(1)/100);  % the series with some randomness

      plot(x,series(counter,2:end),...
      'color',colour) % plot series with a different colour for each class
      
      counter = counter + 1;
    end
end
hold off

% reorder classes randomly:
if reorder_rand == 1
     series_unsorted = series;
    [~, rnd_ind] = sort(rand(1,tot_amt_ser));
    for rows = 1:tot_amt_ser
      series(rows,:) = series_unsorted(rnd_ind(rows),:); 
    end
end

% save series:
if dataset_label == 0;
    save labelled_series.mat series
elseif dataset_label == 1;
    series_train = series;
    save labelled_series_train.mat series_train
elseif dataset_label == 2;
    series_test = series;
    save labelled_series_test.mat series_test
end

disp 'ready'