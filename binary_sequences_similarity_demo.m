%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                        Similarity of Binary Sequences                        %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             doctoral student %
%                                   Music Department - University of Jyväskylä %
%                                                                January, 2016 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested in Octave 4.

% ==============================================================================
% Initialisation:
clc
clear
close all

% ------------------------------------------------------------------------------
% Description:

% A binary sequence is composed by two symbols (e.g. ones and zeroes).

% The following methods measure the similarity of two binary sequences.

% METHOD 1: Pearson's R correlation of the binary sequences convolved with a 
% gaussian kernel (c.f. Hartmann, Lartillot & Toiviainen, 2015).

% METHOD 2: Physical Correlation of the binary sequences convolved with a 
% gaussian kernel (c.f. Bruderer, McKinney & Kohlrausch, 2012).

% METHOD 3: Binary Sequences Similarity

% ------------------------------------------------------------------------------

% Example:
% A participant in an experiment was presented with a music excerpt 10 seconds long.
% The participant was asked to press a button when a change in the music was noticed.
% The button was connected to a machine that recorded the time of pressing the 
% button counting from the beginning of the music excerpt, with a precision of
% 1/10 seconds.
% This means that in total there were 100 potential points where a participant
% could have pressed the button.
% In the first trial the Participant pressed the button at seconds 1.4, 3 and 7.8
% In the second trial the Participant pressed the button at seconds 1, 3.5, 6, 8.4 and 9
% Several methods can be used to measure the similarity of both responses.

% ------------------------------------------------------------------------------
% Enter data and parameters:
% By entering data and parameters below it is possible to investigate the
% implementations of different methods to assess how similar are two binary
% sequences.

sequences_length = 100 % <---- length of the binary sequences

% location of "ones" (all the rest are zeroes) for each binary sequence:
indexes(1) = {[ 14 30    78   ]};
indexes(2) = {[ 10 35 60 84 90]};


% only for methods 1 and 2:

bandwidth = 13; % <----------- Gaussian window (kernel) size 
alpha = 2*pi/bandwidth; % <--- Gaussian window (kernel) alpha factor 

zero_mean = 0; % <------------ make convolved curve normally distributed on/off

% vertical offsets:
v_offset_1 = 0;
v_offset_2 = 0;

% ------------------------------------------------------------------------------
% Report:


disp('--------------------------------------------------')
disp(sprintf('amount of points sequence 1 = %i (%1.2f of the total)',...
   size(cell2mat(indexes(1)),2),(size(cell2mat(indexes(1)),2))/sequences_length ));
disp(sprintf('amount of points sequence 2 = %i (%1.2f of the total)',...
   size(cell2mat(indexes(2)),2),(size(cell2mat(indexes(2)),2))/sequences_length ));
disp('--------------------------------------------------')

% ------------------------------------------------------------------------------
% Make pulse trains by convolving the binary sequences with a gaussian window 
% (methods 1 and 2):

pulse_trains = zeros(2,sequences_length);
gaussian_window = gaussian(bandwidth,alpha);

% Compute Gaussian Kernel (comment/uncomment option a or b).....................

%% a) use Octave's function (from the Signal package):
%pkg load all
%gaussian_window = height * gaussian(bandwidth,alpha);

% b) embedded code:
x = linspace( -(bandwidth-1)/2, (bandwidth-1)/2, bandwidth );
gaussian_window = ( exp( -((alpha*x).^2) / 2 ) ); 

% ..............................................................................

for i = 1:2
   pulse_trains(i,cell2mat(indexes(i))) = 1;
   convoluted_vectors(i,:) = conv(pulse_trains(i,:),gaussian_window);
end

shiftback = round(bandwidth/2);

convoluted_vectors_shifted =  zeros(2,sequences_length);
convoluted_vectors_shifted =  convoluted_vectors(:,shiftback:end);

% ------------------------------------------------------------------------------
% Pre-process data (methods 1 and 2):

% remove mean:
if zero_mean == 1
   convoluted_vectors_shifted(1,:) = convoluted_vectors_shifted(1,:) - mean(convoluted_vectors_shifted(1,:));
   convoluted_vectors_shifted(2,:) = convoluted_vectors_shifted(2,:) - mean(convoluted_vectors_shifted(2,:));
end

% offset vertically:
pulse_trains(1,find(pulse_trains(1,:)==1)) = pulse_trains(1,find(pulse_trains(1,:)==1)) + v_offset_1;
pulse_trains(2,find(pulse_trains(2,:)==1)) = pulse_trains(2,find(pulse_trains(2,:)==1)) + v_offset_2;
convoluted_vectors_shifted(1,:) = convoluted_vectors_shifted(1,:) + v_offset_1;
convoluted_vectors_shifted(2,:) = convoluted_vectors_shifted(2,:) + v_offset_2;

% ==============================================================================
% Compute similarity:
% ------------------------------------------------------------------------------
% METHOD 1: Pearson's R correlation
% Comment/uncomment options 1.a or 1.b

% 1.a - R, external function ...................................................
%Pearsons_R_function = corr(convoluted_vectors_shifted(1,:),convoluted_vectors_shifted(2,:))
%R = Pearsons_R_function;

% 1.b - R, embedded ............................................................
Pearsons_R_embedded = ...
   (...
      ( convoluted_vectors_shifted(1,:) - mean(convoluted_vectors_shifted(1,:)) )...
      *...
      ( convoluted_vectors_shifted(2,:) - mean(convoluted_vectors_shifted(2,:)) )'...
      /... 
      sqrt(...
         sum(( convoluted_vectors_shifted(1,:) - mean(convoluted_vectors_shifted(1,:)) ).^2)...
         *...
         (sum(( convoluted_vectors_shifted(2,:) - mean(convoluted_vectors_shifted(2,:)) ).^2))'...
      )...
   )
R = Pearsons_R_embedded;

% ------------------------------------------------------------------------------
% METHOD 2: Physical Correlation
% Comment/uncomment options 2.a or 2.b

% 2.a - C, external function ...................................................
%Physical_Correlation_C_function = physcorr(convoluted_vectors_shifted(1,:),convoluted_vectors_shifted(2,:))
%C = Physical_Correlation_C_function;

% 2.b - C, embedded ............................................................
Physical_Correlation_C_embedded = (convoluted_vectors_shifted(1,:) * convoluted_vectors_shifted(2,:)')...
   / sqrt( sum(convoluted_vectors_shifted(1,:).^2) * sum(convoluted_vectors_shifted(2,:).^2) )   
C = Physical_Correlation_C_embedded;
   
% ------------------------------------------------------------------------------
% METHOD 3: Binary Sequences Similarity (version 2, year 2016)
% Comment/uncomment options 3.a or 3.b

% 3.a - S, external function ...................................................
% Be sure to use version 2 of year 2016 !!!

%disp('--------------------------------------------------')
%Binary_Sequences_Similarity_S_function = ...
%binseqsi(cell2mat(indexes(1)),cell2mat(indexes(2)),sequences_length,1);
% S = Binary_Sequences_Similarity_S_function;

% 3.b - S, embedded ............................................................
a = cell2mat(indexes(1));
b = cell2mat(indexes(2));
size_a = length(a);
size_b = length(b);
maxsize = max(size_a, size_b);
minsize = min(size_a, size_b);

% make reference matrices for each vector:
id_a = repmat(a',1,size_b);
id_b = repmat(b,size_a,1);

% distance matrix:
for i_1 = 1:size_a %
    for i_2 = 1:size_b
        distmat(i_1,i_2) = -(a(i_1)-b(i_2));
    end
end

absdistmat = abs(distmat);

% Make a logical matrix indicating the minimum value(s) of each column:
for i_1 = 1:size_a
    mincols(i_1,:) = absdistmat(i_1,:) == min(absdistmat);
end

% Make a logical matrix indicating the minimum value(s) of each row:
for i_1 = 1:size_b
    minrows(:,i_1) = absdistmat(:,i_1) == min(absdistmat,[],2);
end

allmins = mincols .* minrows; % intersection
allmins(allmins == 0) = NaN;

% find to which elements of original a and b do these minima correspond:
from_a = allmins.*id_a;
m{1,1} =  (from_a(isfinite(from_a)))';
from_b = allmins.*id_b;
m{2,1} = (from_b(isfinite(from_b)))';

% compute measures:
d = (nanmean(nanmean(allmins(:) .* absdistmat(:)))); % distance
minunel = min(size(unique(m{1,1}),2),size(unique(m{2,1}),2));

c = (1-(d/sequences_length)); % paired elements closeness
f = minunel/maxsize; % fraction of paired elements
Binary_Sequences_Similarity_S_embedded =  c*f
S = Binary_Sequences_Similarity_S_embedded;
   
% ------------------------------------------------------------------------------
% Plot:

pulse_trains_nozeros = pulse_trains;
pulse_trains_nozeros(pulse_trains_nozeros==0) = NaN;

fig = figure('Position',[600,50,600,480]); % set figure position and size

for i = 1:2

   subplot(3,1,i+1)
   plot(pulse_trains_nozeros(i,:),'.','markersize',10,'Color',[0 0 0]+0); % plot pulse train
   hold on
   plot(convoluted_vectors_shifted(i,:),'linewidth',2,'Color',[0 0 0]+0.7); % plot convoluted curve
   set(gca,'xlim',[0,sequences_length]);
   xlabel('samples')
   grid on;
   title(['SEQUENCE ',num2str(i)])
   hold off

end

% put a title to the figure:
title_axes = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
text(0.5,( 1 - 0.28 - 0.00),...
   ['COMPARISON OF SIMILARITY MEASURES FOR BINARY SEQUENCES'],...
   'HorizontalAlignment','center','VerticalAlignment','top')
text(0.5,( 1 - 0.28 - 0.03),...
   ['R = ', num2str(R),...
   ' | C = ', num2str(C),...
   ' | S = ', num2str(S)],...
   'HorizontalAlignment','center','VerticalAlignment','top')

   
% Save the figure:
%saveas(fig,cstrcat(...
%  'COM_SIM_R',num2str(Pearsons_R_builtin),...
%  '_C',num2str(Physical_Correlation_C),...
%  '_S',num2str(Binary_Sequences_Similarity_S),'.png'))
