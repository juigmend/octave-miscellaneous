%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                           Kernel Density Estimation                          %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             doctoral student %
%                                   Music Department - University of Jyväskylä %
%                                                                January, 2016 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested in Octave 4.
% The kernel_density function is part of the econometrics package.

% ==============================================================================
% Initialisation:
clc
clear
close all
pkg load all

% ------------------------------------------------------------------------------
% Description:

% This is a method to estimate density of points in a binary sequence.

% A binary sequence is composed by two symbols (e.g. ones and zeroes).
% Kernel Density Estimation will convolve the points of interest in the 
% sequence (i.e. ones) with a window (often gaussian) and will scale the curve 
% so that the sum of all the points is one.

% The peaks of the curve can be regarded as a representative sequence at
% different density thresholds
% (c.f. Burunat, Alluri, Toiviainen, Numminen & Brattico, 2014).

% ------------------------------------------------------------------------------

% Example:
% Three participants in an experiment were presented with a 10 seconds music excerpt.
% The participants were asked to press a button when a change in the music was noticed.
% The button was connected to a machine that recorded the time of pressing the 
% button counting from the beginning of the music excerpt, with a precision of
% 1/10 seconds.
% This means that in total there were 100 potential points where a participant
% could have pressed the button.
% Participant 1 pressed the button at seconds 1.4, 3 and 7.8
% Participant 2 pressed the button at seconds 1, 3.5, 6, 8.4 and 9 
% Participant 3 pressed the button at seconds 1.2, 4, 5 and 9.6
% Kernel Density Estimation is used to generate data that would best represent
% the combined data of the three participants.

% ------------------------------------------------------------------------------
% Enter data and parameters:
% By entering data and parameters below it is possible to investigate the
% implementation of the method.

sequences_length = 100; % <--- length of raw binary sequences
indexes = [10 12 14 30 35 40 50 60 78 84 90 96]'; % <--- data points

bandwidth = 4; % <------------ Kernel bandwidth, A.K.A window size (KDE)
peak_int = 1; % <------------- Peak threshold (PDE) integer
peak_mult = 10^-2; % <-------- Peak threshold multiplier

% ------------------------------------------------------------------------------
% Make KDE curve:

density_curve = kernel_density([1:sequences_length]',indexes,bandwidth);

% ------------------------------------------------------------------------------
% Extract density peaks:

% Extract peaks by testing for a negative second derivative:
peaksind = find(diff(sign(diff(density_curve))) == -2) + 1;

% Retrieve peaks over a density threshold:
pk = peak_int * peak_mult;
selpeaks = peaksind(find(density_curve(peaksind) >= pk));

% Make a vector of peaks and zeroes (for later plot):
peaks_curve = zeros(1,sequences_length);
for i = 1:size(selpeaks,1)
    peaks_curve(1,(selpeaks(i))) = density_curve(selpeaks(i));
end
peaks_curve(peaks_curve == 0) = NaN;

% Make a peak threshold line:
threshold_line(1,1:sequences_length) = pk;

% ------------------------------------------------------------------------------
% Plot:

fig = figure('Position',[100,10,600,400]); % set figure position and size

subplot(2,1,1)
plot(indexes',1,'.','markersize',10,'Color',[0 0 0]+0); % plot binary sequence (also known as "pulse train")
set(gca,'xlim',[0,sequences_length],'ylim',[0.9,1.1],'yticklabel',[]);
xlabel('samples')
grid on;
title('RAW SEQUENCE')

subplot(2,1,2)
plot([1:sequences_length],density_curve','linewidth',2,'Color',[0 0 0]+0.5); % plot estimated density curve
set(gca,'xlim',[0,sequences_length]);
xlabel('samples')
ylabel('density')
grid on;
hold on;
plot([1:sequences_length],peaks_curve,'.','markersize',10,'Color',[0 0 0]+0); % plot density peaks
plot(threshold_line,'-','LineWidth',4,'Color',[0 0 0]+0.8) % plot peak threshold line
title(strcat('DENSITY ESTIMATION | KDE = ',num2str(bandwidth),', PDE = ',num2str(pk)))

% Save the figure:
%saveas(fig,strcat('KDE_TEST_8ve-',num2str(bandwidth),'-',num2str(peak_int),'.png'))
%saveas(fig,strcat('KDE_TEST_8ve-',num2str(bandwidth),'-',num2str(peak_int),'.pdf'))