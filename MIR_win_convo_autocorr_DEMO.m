%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                                                              %
%                          MUSIC INFORMATION RETRIEVAL                         %
%                   WINDOWING, CONVOLUTION AND AUTOCORRELATION                 %
%                                                                              %
%                                 September, 2017                              %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested with Matlab R2015a

% ==============================================================================
% Instructions:

%    Run the cells one by one, following the comments. The beginning of a cell is 
%    indicated by a double percentage sign (%%) to the left of a line.
%    A cell may be run by depressing keys 'cmd' and 'enter' at the same time.

%...............................................................................
% Contents:

% line  42: 1) WINDOWING
% line 174: 2) CONVOLUTION 
% line 221: 3) CORRELATION, CROSS-CORRELATION AND AUTOCORRELATION
% line 291: 4) MUSIC INFORMATION RETRIEVAL EXAMPLE: Simple algorithm for pitch contour extraction
% line 436: BONUS 1) Estimation of fundamental frequency by extracting harmonics using Fast Fourier Transform
% line 488: BONUS 2) Segmentation by convolving a checkerboard kernel upon the diagonal of a distance matrix

%% ==============================================================================
% Initialisation:

clc
clear
close all

%% =============================================================================
% 1) WINDOWING

% "Windowing" refers to a process in which computation is performed upon a section
% (a "window") of a signal. This process can be repeated at several and 
% equally-spaced points along the signal.

% The windows are also referred to as "frames".

% Make a short test signal (a wavelet):
wavelet_1 = [0 1 2 3 4 5 4 3 2 2 3 4 3 2 1 0 ];
plot(wavelet_1)

%%
% Make the wavelet negative:

wavelet_2 = wavelet_1 * -1;
plot(wavelet_2)

%%
% Put the positive and negative wavelets side by side:

wavelet_3 = [wavelet_1, wavelet_2];
plot(wavelet_3)

%%
% The Homer Simpson method of making a signal (copy and paste the wavelet many times):

test_signal_1 = [wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3, wavelet_3];
plot(test_signal_1)

%%
% The Matlab method (use the function 'repmat'):

test_signal_2 = repmat(wavelet_3,1,20);
plot(test_signal_2)

%%
% Extract a window from the signal:

single_window = test_signal_2(12:30);
plot(single_window)

%%
% Make a window so that its boundaries can be specified with variables:

window_width = 18;
window_start = 12;
window_end = window_start + window_width - 1 ;
single_window = test_signal_2(window_start:window_end);
plot(single_window)

%%
% Make a sliding window with hop = 1 sample:

window_width = 2; 
amount_windows = length(test_signal_2) - window_width - 1 ;

processed = zeros(1,amount_windows);

for i = 1:amount_windows
    window_start = i;
    window_end = window_start + window_width;
    processed(1,i) = mean(test_signal_2(window_start:window_end)); % process is mean (moving average)
%     processed(1,i) = rms(test_signal_2(window_start:window_end)); % process is root-mean-square (a measure of power)
%     processed(1,i) = std(test_signal_2(window_start:window_end)); % process is standard deviation (difference from mean)
end

plot(processed)

% Observe the results at different values for window width and process.

%%
% Make a sliding window with assignable hop:

window_width = 8; 
hop = 2;
amount_windows = round( (length(test_signal_2) - window_width) / hop );

processed = zeros(1,amount_windows);

for i = 1:amount_windows
    window_start = ((i - 1) * hop) + 1; % incorporate the hop to the indexing
    window_end = window_start + window_width;
    processed(1,i) = mean(test_signal_2(window_start:window_end)); % process is mean (moving average)
%     processed(1,i) = rms(test_signal_2(window_start:window_end)); % process is root-mean-square (a measure of power)
%     processed(1,i) = std(test_signal_2(window_start:window_end)); % process is standard deviation (difference from mean)
end

plot(processed)

% Observe the results at different values for window width, hop and process.

%%
% Upsample the test signal:

resampling_rate = 20;
test_signal_3 = zeros(1,length(test_signal_2) * resampling_rate );

counter = 1;
for i = 1:length(test_signal_2)
    for i_1 = 1:resampling_rate
        test_signal_3(counter) = test_signal_2(i);
        counter = counter + 1;
    end
end

plot(test_signal_3)

%%
% Run sliding window process upon upsampled test signal:

tic % start stopwatch

window_width = 100; 
hop = 10;
amount_windows = round( (length(test_signal_3) - window_width) / hop );

processed = zeros(1,amount_windows);

for i = 1:amount_windows
    window_start = ((i - 1) * hop) + 1; % incorporate the hop to the indexing
    window_end = window_start + window_width - 1 ;
    processed(1,i) = mean(test_signal_3(window_start:window_end)); % process is mean (moving average)
end

toc % stop stopwatch

plot(processed)

% Observe the difference in computing time using different values for window width and hop.

%% =============================================================================
% 2) CONVOLUTION

% Convolution can be thought of as a special kind of windowed process, in which
% the windowed section of the signal is dot-multiplied with the reversed version
% of another signal that is the same size of the window. This second signal is
% often much shorter than the first one. In Computer Science jargon the process
% performed over each window along the signal, is called "kernel".
% Also usually the second (shorter) signal itself is called "kernel".

% Make kernel:

kernel_signal = [0 1 0.7 0.5 0.2 0]';
% kernel_signal = gausswin(20);
plot(kernel_signal)

%%
% Flip (reverse) the kernel:

flipped_kernel_signal = flipud(kernel_signal);
plot(flipped_kernel_signal)

%%
% Run sliding window with hop = 1 sample, where process is dot product:

signal_to_process = test_signal_2;

window_width = length(kernel_signal); 
hop = 1; % for convolution hop should be 1 
amount_windows = round( (length(signal_to_process) - window_width) / hop );

processed = zeros(1,amount_windows);

for i = 1:amount_windows
    window_start = ((i - 1) * hop) + 1; % incorporate the hop to the indexing
    window_end = window_start + window_width - 1 ;
    processed(1,i) = signal_to_process(window_start:window_end) * flipped_kernel_signal; % process is dot product
end

plot(processed)

% The previous examples show convolution of two one-dimensional signals.
% Convolution can also be performed in more than one dimension.
% For example, two-dimensional convolution it is widely used to filter images.
% Convolution in two dimensions takes the dot product of a matrix and another
% matrix, the latter usually being smaller.

%% =============================================================================
% 3) CORRELATION, CROSS-CORRELATION AND AUTOCORRELATION

% Correlation is a measure of the relationship between two variables.
% It can be considered as a measure of similarity.
% Intuitively, there are three cases:

% 1) Full correlation: Two variables A and B change in the same proportion.
%    This can be expressed as corr(A,B) = 1
%    Hence, always corr(A,A) = 1

% 2) Negative correlation: Two variables A and B change in the same inverse 
%    proportion.
%    This can be expressed as corr(A,B) = -1
%    Hence, always corr(A,-A) = -1

% 3) No correlation: When either variable A or B change, the other variable
%    doesn't change. This is undefined (there is no value for it).

% A common measurement of correlation is Pearson's R, which consist of a rather 
% lengthy formula:

A = [1 2 3 4];
B = [1 2 3 4];

R = ( ( B - mean(B) ) *  ( A - mean(A) )' ) / sqrt( sum( (A - mean(A)).^2 ) * sum( (B - mean(B)).^2 ));
R = corr(A',B');

%%

% Cross-correlation is a much simpler operation. The concept is still of 
% similarity between two signals. It is obtained by computing the dot product of
% the signals at different displacements (lags).
% The resulting value will be the highest at the displacement point (lag) where 
% the signals are most similar, and viceversa.
% Computationally, it can be seen as a windowed process, actually the same as
% convolution except that the second signal to be multiplied is not reversed.
% Also, although convolution and cross-correlation are mathematically very similar,
% conceptually they are quite different. Convolution is a process in which the 
% longer signal acquires characteristics of the shorter signal. This is, for 
% example, how convolution reverberation and filters (such as audio equalizers) 
% work. On the other hand, cross-correlation is a measure of similarity between 
% two signals, at different displacements.

% Auto-correlation is simply cross-correlation of a signal with itself.

% Run sliding window with hop = 1 sample, where process is dot product:

signal_to_process_1 = test_signal_3;
signal_to_process_2 = test_signal_3;

window_width = length(signal_to_process_2); 
hop = 1;
amount_windows = length(signal_to_process_1)*2 - 1;

processed = zeros(1,amount_windows);

% Add zeros at the beginning and end of one of the signals, to slide the second
% signal completely:
signal_to_process_1_padded = zeros(1,length(signal_to_process_1)*3);
signal_to_process_1_padded(window_width:window_width+length(signal_to_process_1)-1) = signal_to_process_1;

for i = 1:amount_windows
    window_start = ((i - 1) * hop) + 1; % incorporate the hop to the indexing
    window_end = window_start + window_width - 1 ;
    processed(1,i) = signal_to_process_1_padded(window_start:window_end) * signal_to_process_2'; % process is dot product
end

plot(processed)

%% =============================================================================
% 4) MUSIC INFORMATION RETRIEVAL EXAMPLE: 
%    Simple algorithm for pitch contour extraction

example_files_path = '/xxxx/xxxx'; % <--------------- Path where the example files are
addpath(example_files_path)

%%
% Load audio file :
audio_filename = 'xxxx'; % e.g., 'man_sings_note_C.wav'
[audio_raw, audio_sf] = audioread(audio_filename);

% Play raw audio file:
sound(audio_raw, audio_sf)

%%
% Extract and play window:

window_width = audio_sf * 0.05;
window_start = audio_sf * 0.3;
window_end = window_start + window_width;

audio_window = audio_raw(window_start:window_end);
sound(audio_window, audio_sf)

%%
% Estimation of fundamental frequency by autocorrelation:

amount_windows = length(audio_window)*2;
processed_full = zeros(amount_windows-1,1);
audio_window_padded = zeros(length(audio_window)*3,1);
audio_window_padded(window_width:window_width+length(audio_window)-1) = audio_window;
for i = 1:amount_windows-1
    window_start = i;
    window_end = window_start + window_width;
    processed_full(i) = audio_window_padded(window_start:window_end)' * audio_window;
end
processed = processed_full(floor(window_width/2)+1:floor(window_width*1.5)); 

figure; plot(processed)

% Note that the result is symmetric.

%%
% Intuitively try to extract highest peaks:

maximum_index = find(processed == max(processed)); % maximum peak
% window_to_examine = processed(maximum_index-1:maximum_index+201);
window_to_examine = processed(maximum_index-1:end);

figure; 
% plot(processed)
plot(window_to_examine )
hold on
grid on
% plot( diff(window_to_examine) ,'linewidth',3 )
% plot( sign(diff(window_to_examine)) ,'linewidth',3 )
% plot( diff(sign(diff(window_to_examine))) ,'linewidth',3 )
plot( diff(sign(diff(window_to_examine))) < 0,'linewidth',3 )

%%
% Compute the period (T, for time) between the highest peaks and extract frequency (f):
% Remember that T = 1/f ; hence f = 1/T

selected_peaks_bool = diff(sign(diff(window_to_examine))) < 0;
values_selpeaks = selected_peaks_bool .* window_to_examine(1:end-2);
second_biggest_peak = max(values_selpeaks(2:end));
T = find(window_to_examine == second_biggest_peak); % period in samples
freq = 1/T;             % frequency in samples 
freq = freq * audio_sf; % frequency in Hz

%%
% Run windowed frequency estimation along an audio signal.
% This will result in a 'pitch contour'
% This process invoves a windowed process inside another windowed process.
% The outer loop will slide the windowed pitch estimation process along the audio file.
% The inner loop will slide a dot product along the window with the window itself and
% will compute the frequency.

% Load audio file :
audio_filename = 'xxxx'; % e.g., 'child_sings_scale_C_major.wav'
[audio_raw, audio_sf] = audioread(audio_filename);

% Play raw audio file:
sound(audio_raw, audio_sf)

%%
clc
close all
tic

margin = round(audio_sf / 10);
audio_raw = audio_raw(margin:end - margin);

window_width = audio_sf * 0.01; 
window_width = round(window_width);
hop = round(audio_sf / 40);
outer_amount_windows = round( (length(audio_raw) - window_width) / hop );
frequency_contour = zeros(1,outer_amount_windows);

inner_amount_windows = window_width;
inner_autocorr = zeros(inner_amount_windows,1);
audio_window_padded = zeros(window_width * 2,1);

for i = 1:outer_amount_windows
    
    outer_window_start = ((i - 1) * hop) + 1; % incorporate the hop to the indexing
    outer_window_end = outer_window_start + window_width - 1 ;
    audio_window = audio_raw(outer_window_start:outer_window_end);
    
%     pause(0.1)
%     plot(outer_window)

    %...........................................................................
    % estimate fundamental frequency by autocorrelation ('inner' process):
   
    audio_window_padded(1:window_width) = audio_window;
    
%     pause(0.1)
%     plot(audio_window_padded)
    
    for i_1 = 1:inner_amount_windows
        inner_window_start = i_1;
        inner_window_end = inner_window_start + window_width-1;
        inner_autocorr(i_1) = audio_window_padded(inner_window_start:inner_window_end)' * audio_window;
    end
    
%     pause(0.2)     
%     plot(inner_autocorr)

    %...........................................................................
    
    selected_peaks_bool = diff(sign(diff(inner_autocorr))) < 0;
    values_selpeaks = selected_peaks_bool .* inner_autocorr(1:end-2);
    second_biggest_peak = max(values_selpeaks(2:end));
    T = find(inner_autocorr == second_biggest_peak);        % period in samples
    frequency_contour(i) = 1/T;                             % frequency in samples
    frequency_contour(i) = frequency_contour(i) * audio_sf; % frequency in Hz
    
end

toc

plot(frequency_contour)

%% =============================================================================
% BONUS 1)  Estimation of fundamental frequency by extracting harmonics using Fast Fourier Transform

% Load audio file :
audio_filename = 'xxxx'; % e.g., 'man_sings_note_C.wav'
[audio_raw, audio_sf] = audioread(audio_filename);

% Play raw audio file:
sound(audio_raw, audio_sf)

%%
% Compute FFT upon a window:

audio_window = audio_raw(audio_sf/100:audio_sf*2);
audio_window_fft = fft(audio_window,audio_sf);
audio_window_fft_abs = abs(audio_window_fft);
plot(audio_window_fft_abs)

%%
% Compute spectrum from FFT:
% (just take the half of FFT and rescale so maximum is 1)

fft_length = length(audio_window_fft_abs);
audio_window_fft_abs_spectrum = audio_window_fft_abs( 1 : floor(fft_length/2) + 1 );
audio_window_fft_abs_spectrum = audio_window_fft_abs_spectrum / max(audio_window_fft_abs_spectrum); % normalise it, don't criticise it
plot(audio_window_fft_abs_spectrum)

%%
% Set a frequency scale if needed:

frequency_scale = audio_sf * ([0 : floor(fft_length/2 )] / fft_length ) ;
plot(frequency_scale,audio_window_fft_abs_spectrum); 

%%
% Visualise it better in a logarithmic scale:

figure; semilogx(frequency_scale,audio_window_fft_abs_spectrum); 

%%
% Estimation of fundamental pitch by extracting harmonics:

h_threshold = 0.38;
harmonics = (audio_window_fft_abs_spectrum > h_threshold);
semilogx(frequency_scale,audio_window_fft_abs_spectrum); 
hold on
semilogx(frequency_scale,harmonics,'r');
f0 = frequency_scale(find(harmonics==1,1));

% This pitch estimation process can also be applied as a sliding window. 
% The FFT alone applied as a windowed process and then visualised is called 
% 'spectrogram'.

%% =============================================================================
% BONUS 2) Segmentation by convolving a checkerboard kernel upon the diagonal of a distance matrix

% B2.1) Self-Similarity Matrix 

% A self-similarity matrix is a distance matrix. It contains the distance 
% (e.g., the absolute difference) between every sample of a finite signal with
% itself and every other samples.
% Example:

my_signal = [1 2 3 4];

my_distance_matrix_1 = abs( [ 1-1, 1-2, 1-3, 1-4 ;
                              2-1, 2-2, 2-3, 2-4 ; 
                              3-1, 3-2, 3-3, 3-4 ;
                              4-1, 4-2, 4-3, 4-4 ] )

close all
imagesc(my_distance_matrix_1)

%%                          
% This can be quickly computed with the Matlab command 'pdist2'                       

my_distance_matrix_1 = pdist2(my_signal',my_signal')
imagesc(my_distance_matrix_1)

%%
% Compute a distance matrix of a signal with itself (a 'self-similarity matrix'):

my_distance_matrix_2 = pdist2(test_signal_2',test_signal_2');
imagesc(my_distance_matrix_2)

%%
% Compute the self-similarity matrix for the frequency contour of the MIR example:

clc
close all

figure; plot(frequency_contour);

frequency_contour_dm = pdist2(frequency_contour',frequency_contour');
figure; imagesc(frequency_contour_dm)

% %%
% To make the notes more clear we are smooth the frequency contour
% by convolving it with a filter window.
% Now that we know how convolution works, we can use the built-in function 'conv' :)

clc
close all

figure; plot(frequency_contour);
title('frequency contour')

% Try different windows and windows lengths:
window_length = 4;
% filter_window = gausswin(window_length);
% filter_window = hann(window_length);
filter_window = hamming(window_length);
figure; plot(filter_window);
title('Hann window')

smooth_frequency_contour = conv(frequency_contour',filter_window,'same');

figure; plot(smooth_frequency_contour);
title('smooth frequency contour')

smooth_frequency_contour_dm = pdist2(smooth_frequency_contour,smooth_frequency_contour);
figure; imagesc(smooth_frequency_contour_dm)

%%
% B2.2) Novelty Score

% A novelty score can be extracted from a self-similarity simply by convolving a 
% checkerboard kernel along the diagonal of the self-similarity matrix.

% Make a 2x2 checkerboard kernel:

checkerboard_kernel_1 = [ -1,  1 ;
                           1, -1 ]
                       
imagesc(checkerboard_kernel_1)

%%
% Make a checkerboard kernel of variable size:

ck_size = 20;
checkerboard_kernel_2 = [ -ones(ck_size/2),  ones(ck_size/2) ;
                           ones(ck_size/2), -ones(ck_size/2) ];                       
imagesc(checkerboard_kernel_2)

%%
% To convolve the checkerboard kernel along the diagonal of the self-similarity 
% matrix the easy way is to just convolve the whole matrix and then extract the
% diagonal. This is easy to code but computationally inefficient.

% Convolve the kernel and the self-similarity matrix using the built-in
% function 'conv2':

ssm_convolved = conv2(smooth_frequency_contour_dm,checkerboard_kernel_2,'same');
imagesc(ssm_convolved)

%%
% Extract the diagonal of the convolved matrix:

novelty_score = diag(ssm_convolved);
figure; plot(novelty_score)

%%
% The peaks of the novelty score are the boundaries of each note, so next we find
% where those peaks are. Now that we know how peak finding works we can simply 
% use the built-in function 'findpeaks' :D
clc
close all

[peaks_values, peaks_indexes] = findpeaks(novelty_score);

peaks_bool = zeros(1,length(smooth_frequency_contour));
peaks_bool(peaks_indexes) = 1;
peaks_to_plot = peaks_bool * max(smooth_frequency_contour);

plot(smooth_frequency_contour)
hold on
plot(peaks_to_plot)

%%
% B2.3) Extract average frequency of each segment

amount_of_segments = length(peaks_indexes);
boundaries_indexes = [1; peaks_indexes]; % we add a one at the beginning to indicate the first boundary
mean_freqs = zeros(amount_of_segments,1);

for i = 1:amount_of_segments
    
   start_of_segment = boundaries_indexes(i);
     end_of_segment = boundaries_indexes(i+1);
     
   mean_freqs(i) = mean( frequency_contour( start_of_segment : end_of_segment) );
   
end

mean_freqs









