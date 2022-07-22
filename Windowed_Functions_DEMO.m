%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %
% Windowed Functions DEMO %
%                         %
%               July,2017 %
%    Juan Ignacio Mendoza %
% University of Jyväskylä %
%                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tested in GNU-Octave 4

clc
clear
close all

line_width = 1;
screen_div = 1.5;
screen_offset = [260, 110];
screen_size = get(0,'screensize');
figures_position(3:4) = screen_size(3:4) ./ screen_div ;
figures_position(1:2) = screen_offset ;
%figure('position', figures_position)

% ==============================================================================
%                   Windowed functions in one dimenson
%% -----------------------------------------------------------------------------
% MAKE ARTIFICIAL TIME SERIES

art_tseries_length = 1200;
art_tseries_max = 1; 
art_tseries_min = -1;
amplitude = (art_tseries_max + abs(art_tseries_min))/2;
period = 400; 
frequency = 1/period ;
offset = period/4;
%period = 1/frequency;

% one in the middle of zeroes:
art_tseries.point = zeros(1,art_tseries_length);
art_tseries.point(round(art_tseries_length/2) + fix(offset)) = 1;

% diversely-sized clusters of ones:
art_tseries.oneclusters = zeros(1,art_tseries_length);
amount_clusters = fix(art_tseries_length/period);
max_elements_cluster = 5;
amt_elements_clusters = round(rand(1,amount_clusters) * max_elements_cluster);
max_size_cluster = period;
sizes_clusters = round(rand(1,amount_clusters) * max_size_cluster);
position_clusters = round((rand(1,amount_clusters) * art_tseries_length + fix(offset)));
for i = 1:amount_clusters
  clusters{i} = round( rand(1,amt_elements_clusters(i)) * sizes_clusters(i) ) + position_clusters(i);
  art_tseries.oneclusters(clusters{i}) = 1;
end

% [1 1 1 1 1 0 0 0 0 0] repeat this pattern:
square_max = repmat(art_tseries_max,1,period/2);
square_min = repmat(art_tseries_min,1,period/2);
art_tseries.square = repmat([square_max,square_min],1,fix(art_tseries_length/period)+2);
art_tseries.square = art_tseries.square(offset+1:art_tseries_length+offset+1);

% sinusoid:
art_tseries.sine_art_tseries = amplitude * sin( 2 * pi * frequency * [0:1:(art_tseries_length-1)] + offset*2*pi/period);

% frequency modulated sinusoid:
fm_period = period/4; 
fm_offset = offset; 
fm_frequency = 1/fm_period ;
fm_strength = 2;
fm_message = fm_strength * sin( 2 * pi * fm_frequency * [0:1:(art_tseries_length-1)] + fm_offset*2*pi/period);
art_tseries.fm = amplitude * sin( 2 * pi * frequency * [0:1:(art_tseries_length-1)] + fm_message);

% amplitude frequency modulated sinusoid:
amfm_factor = 2;
am_period = fm_period*amfm_factor; 
am_offset = offset;
am_frequency = 1/am_period ;
am_strength = 1;
am_message = am_strength * sin( 2 * pi * am_frequency * [0:1:(art_tseries_length-1)] + am_offset*2*pi/period);
art_tseries.amfm = (amplitude + am_message) .* art_tseries.fm;

% add a little noise:
noise_amplitude = amplitude;
noise = rand(1,art_tseries_length);
noise = noise / max(noise');
art_tseries.amfm = art_tseries.amfm + noise_amplitude * noise;

%art_tseries.test_art_tseries = art_tseries.point;
%art_tseries.test_art_tseries = art_tseries.square;
%art_tseries.test_art_tseries = art_tseries.sine_art_tseries;
%art_tseries.test_art_tseries = art_tseries.fm;
%art_tseries.test_art_tseries = art_tseries.amfm;
art_tseries.test_art_tseries = art_tseries.amfm;

% plot artificial time series:
art_tseries_names = fieldnames(art_tseries);
amount_ats_names = length(art_tseries_names);
figure('position', figures_position)
for i = 1:amount_ats_names
  subplot(amount_ats_names,1,i)
  plot(art_tseries.(art_tseries_names{i}),'linewidth',line_width)
  xlim([1,art_tseries_length]);
  ylim([min(art_tseries.(art_tseries_names{i})), max(art_tseries.(art_tseries_names{i}))])
  set(gca,'ytick',([min(art_tseries.(art_tseries_names{i})), 0, max(art_tseries.(art_tseries_names{i}))]))
  title(art_tseries_names{i},'interpreter','none')
end

%% -----------------------------------------------------------------------------
% MAKE A FUNCTION WINDOW
% Window may be also called frame

window_size = art_tseries_length;
%window_size = 10;
window_start = 1;
window_end = window_start + window_size - 1;
function_window = zeros(1,window_size);
function_window = art_tseries.test_art_tseries(1,window_start:window_end);
function_output = mean(function_window); % <--------------- the mean function :(

% plot window;
figure('position', figures_position)
subplot(2,1,1)
plot(function_window,'linewidth',line_width)
ylim([min(art_tseries.(art_tseries_names{i})), max(art_tseries.(art_tseries_names{i}))])
set(gca,'ytick',([min(art_tseries.(art_tseries_names{i})), 0, max(art_tseries.(art_tseries_names{i}))]))
title('windowed raw art_tseries')
subplot(2,1,2)
plot(function_output,'linewidth',line_width)
title('function output')

%% -----------------------------------------------------------------------------
% MAKE A SLIDING FUNCTION WINDOW

%window_size = art_tseries_length;
%window_size = art_tseries_length/10;
window_size = 20;

% add zero-pad margins to a time series:
extended_art_tseries.test_art_tseries_zeropad = zeros(1,art_tseries_length + window_size);
window_margin = window_size/2; 
extended_art_tseries.test_art_tseries_zeropad( window_margin + 1 : (art_tseries_length + window_margin)) = art_tseries.test_art_tseries;

% add border-element padding:
border_elements(1) = art_tseries.test_art_tseries(1);
border_elements(2) = art_tseries.test_art_tseries(end);
extended_art_tseries.test_art_tseries_borderpad = extended_art_tseries.test_art_tseries_zeropad;
extended_art_tseries.test_art_tseries_borderpad(1:window_margin) = border_elements(1);
extended_art_tseries.test_art_tseries_borderpad(window_margin+art_tseries_length:end) = border_elements(2);

% add inner-margin-average padding:
extended_art_tseries.test_art_tseries_meanpad = extended_art_tseries.test_art_tseries_zeropad; % init the mean pad :(
extended_art_tseries.test_art_tseries_meanpad(1:window_margin) = mean(art_tseries.test_art_tseries(1:window_margin));
extended_art_tseries.test_art_tseries_meanpad(window_margin+art_tseries_length:end) = mean(art_tseries.test_art_tseries(art_tseries_length-window_margin+1:end));

extended_art_tseries.test_art_tseries_padded = extended_art_tseries.test_art_tseries_meanpad; % select padding

% plot extended time series:
extended_art_tseries_names = fieldnames(extended_art_tseries);
amount_extended_art_tseries = length(extended_art_tseries_names);
figure('position', figures_position)
for i = 1:amount_extended_art_tseries
  subplot(amount_extended_art_tseries,1,i)
  plot(extended_art_tseries.(extended_art_tseries_names{i}),'linewidth',line_width)
  xlim([1,length(extended_art_tseries.(extended_art_tseries_names{i}))]);
  ylim([min(extended_art_tseries.(extended_art_tseries_names{i})), max(extended_art_tseries.(extended_art_tseries_names{i}))])
  set(gca,'ytick',([min(extended_art_tseries.(extended_art_tseries_names{i})), 0, max(extended_art_tseries.(extended_art_tseries_names{i}))]))  
  title(extended_art_tseries_names{i},'interpreter','none')
end

% compute sliding window:
window_functions.test_art_tseries = art_tseries.test_art_tseries; 
window_functions.window_function_padded = zeros(size(extended_art_tseries.test_art_tseries_padded)); % init
for i = 1 : art_tseries_length
  window_start = i;
  window_end = window_start + window_size - 1;
  function_window = zeros(1,window_size);
  function_window = extended_art_tseries.test_art_tseries_padded(1,window_start:window_end);
  function_output = mean(function_window); % <------------- the mean function :(
  window_functions.window_function_padded(1,i + window_margin) = function_output;
end

% output of same length (get rid of outer margins):
window_functions.window_function_same = window_functions.window_function_padded(window_margin +1 : end - window_margin); 

% output of valid length (get rid of outer and inner margins):
window_functions.window_function_valid = window_functions.window_function_same(window_margin +1 : end - window_margin); 

% plot applied window functions in different shapes (padded, same, valid):
% padded a.k.a. 'full'
window_functions_names = fieldnames(window_functions);
amount_wf_names = length(window_functions_names);
figure('position', figures_position)
for i = 1:amount_wf_names
  subplot(amount_wf_names,1,i)
  plot(window_functions.(window_functions_names{i}),'linewidth',line_width)
  xlim([1,length(window_functions.(window_functions_names{i}))])
  ylim([min(window_functions.(window_functions_names{i})), max(window_functions.(window_functions_names{i}))])
  set(gca,'ytick',([min(window_functions.(window_functions_names{i})), 0, max(window_functions.(window_functions_names{i}))]))  
  title(window_functions_names{i},'interpreter','none')
end

%% -----------------------------------------------------------------------------
% MAKE DIFFERENT KERNELS

kernel_width = window_size;

% square window:
kernels.square = ones(1,kernel_width);

% triangle window:
kernels.triangle = kernels.square;
kernels.triangle(1:fix(kernel_width/2)) = 0 : 1 / ( fix(kernel_width/2) - 1 ) : 1 ;
kernels.triangle(fix(kernel_width/2)+1:end) = flip(kernels.triangle(1:fix(kernel_width/2))) ;

% trapezoidal window:
kernels.trapezoid = kernels.square;
trapezoid_tilt = kernel_width/10;
kernels.trapezoid(1:fix(trapezoid_tilt)) = 0 : 1 / ( fix(trapezoid_tilt) - 1 ) : 1 ;
kernels.trapezoid(end - fix(trapezoid_tilt)+1:end) = flip(kernels.trapezoid(1:fix(trapezoid_tilt))) ;

% Hamming window:
kernels.hamming = [0:kernel_width-1];
hamming_alpha = 25/46;
hamming_beta = 1 - hamming_alpha;
kernels.hamming = hamming_alpha - hamming_beta * cos(2*pi*kernels.hamming/(kernel_width-1));

% Gaussian unit area:
gauss_alpha = 2*pi/kernel_width;
gauss_height = (sqrt(2*pi))/kernel_width;
x = linspace( -(kernel_width-1)/2, (kernel_width-1)/2, kernel_width );
kernels.gauss = gauss_height * ( exp( -((gauss_alpha *x).^2) / 2 ) );

kernels.kernel = kernels.gauss; % select kernel

% plot 1D-kernels:
kernel_names = fieldnames(kernels);
amount_kernel_names = length(kernel_names);
figure('position', figures_position)
for i = 1:amount_kernel_names
  subplot(amount_kernel_names,1,i)
  plot(kernels.(kernel_names{i}),'linewidth',line_width)
  xlim([1,kernel_width])
  set(gca,'ytick',([min(kernels.(kernel_names{i})), 0, max(kernels.(kernel_names{i}))]))
  title(kernel_names{i},'interpreter','none')
end

%% -----------------------------------------------------------------------------
% CONVOLUTION
% Multiply piecewise (dot product) and sum the kernel and each window

convolutions.test_art_tseries = art_tseries.test_art_tseries; 
convolutions.convolution_padded = zeros(size(extended_art_tseries.test_art_tseries_padded)); % init

for i = 1 : art_tseries_length
  window_start = i;
  window_end = window_start + kernel_width - 1;
  function_window = zeros(1,kernel_width);
  function_window = extended_art_tseries.test_art_tseries_padded(1,window_start:window_end);
  function_output = sum(function_window.*kernels.kernel); % <--- the convolution function :)
  convolutions.convolution_padded(1,i + window_margin) = function_output;
end

% output of same length (get rid of outer margins):
convolutions.convolution_same = convolutions.convolution_padded(window_margin +1 : end - window_margin); 

% output of valid length (get rid of outer and inner margins):
convolutions.convolution_valid = convolutions.convolution_same(window_margin +1 : end - window_margin); 

% plot convolution in different shapes (padded, same, valid):
convolution_names = fieldnames(convolutions);
amount_wf_names = length(convolution_names);
figure('position', figures_position)
for i = 1:amount_wf_names
  subplot(amount_wf_names,1,i)
  plot(convolutions.(convolution_names{i}),'linewidth',line_width)
  xlim([1,length(convolutions.(convolution_names{i}))])
  ylim([min(convolutions.(convolution_names{i})), max(convolutions.(convolution_names{i}))])
  set(gca,'ytick',([min(convolutions.(convolution_names{i})), 0, max(convolutions.(convolution_names{i}))]))  
  title(convolution_names{i},'interpreter','none')
end

%% -----------------------------------------------------------------------------
% MAKE A SLIDING FUNCTION WINDOW WITH HOP

hop = fix(window_size*2);
%hop = 1;
amt_hopwindows = fix(art_tseries_length/hop);
hopwindows_margin = fix(window_margin/hop);
hopslides.test_art_tseries = art_tseries.test_art_tseries; 
hopslides.hopslide_padded = zeros(1,amt_hopwindows+window_size/hop); % init

for i = [1:amt_hopwindows]
  window_start = i*hop - hop + 1;
  window_end = window_start + window_size - 1;
  function_window = zeros(1,kernel_width);
  function_window = extended_art_tseries.test_art_tseries_padded(1,window_start:window_end);

%  function_output = mean(function_window)
  function_output = sum(function_window.*kernels.triangle);
%  function_output = sum(function_window.*kernels.gauss); 
 
  hopslides.hopslide_padded(1,i+hopwindows_margin) = function_output;
end

% output of same length (get rid of outer margins):
hopslides.hopslide_same = hopslides.hopslide_padded( hopwindows_margin +1 : end - hopwindows_margin ); 

% output of valid length (get rid of outer and inner margins):
hopslides.hopslide_valid = hopslides.hopslide_same( hopwindows_margin +1 : end - hopwindows_margin ); 

% plot processed window functions with hop:
hopslides_names = fieldnames(hopslides);
amount_wf_names = length(hopslides_names);
figure('position', figures_position)
for i = 1:amount_wf_names
  subplot(amount_wf_names,1,i)
  plot(hopslides.(hopslides_names{i}),'linewidth',line_width)
  xlim([1,length(hopslides.(hopslides_names{i}))])
  ylim([min(hopslides.(hopslides_names{i})), max(hopslides.(hopslides_names{i}))])
  set(gca,'ytick',([min(hopslides.(hopslides_names{i})), 0, max(hopslides.(hopslides_names{i}))]))  
  title(hopslides_names{i},'interpreter','none')
end

%% -----------------------------------------------------------------------------
% APPLY WINDOWED FUNTION TO AN AUDIO TIME SERIES

addpath(genpath('D:\Octave code\'));
audio_filename = 'D:\Octave code\Mulla_Sanat_On_puolikertosae.wav';
[audiosignals.raw, audio_sample_frequency] = audioread(audio_filename);
audiosignals.raw = ( audiosignals.raw(:,1) + audiosignals.raw(:,2) ) / 2; % mono mix

% Play original:
sound(audiosignals.raw, audio_sample_frequency)

% Implement computation of windowed function as a function (m.file):
filter_frequency = audio_sample_frequency/60; % Hz
audio_window_size = audio_sample_frequency/filter_frequency ;
%audio_hop = window_size;
audio_hop = 1;
tic
audiosignals.trimmed = audiosignals.raw(1:audio_sample_frequency*2,:);
audiosignals.processed = window_function(audiosignals.trimmed', audio_sample_frequency, 'gauss_kernel', audio_window_size, audio_hop);

% measure computation time (try different window sizes and hops):
display(sprintf('AUDIO PROCESSING TIME = %f s.', toc));
audiosignals.processed = audiosignals.processed';
%whos audiosignals.processed 

% plot audio signals:
audiosignals_names = fieldnames(audiosignals);
amount_as_names = length(audiosignals_names);
figure('position', figures_position)
for i = 1:amount_as_names
  subplot(amount_as_names,1,i)
  plot(audiosignals.(audiosignals_names{i}),'linewidth',line_width)
  xlim([1,length(audiosignals.(audiosignals_names{i}))])
  ylim([min(audiosignals.(audiosignals_names{i})), max(audiosignals.(audiosignals_names{i}))])
  set(gca,'ytick',([min(audiosignals.(audiosignals_names{i})), 0, max(audiosignals.(audiosignals_names{i}))]))  
  title(audiosignals_names{i},'interpreter','none')
end 

% Play processed: 
sound(audiosignals.trimmed, audio_sample_frequency )
sound(audiosignals.processed, audio_sample_frequency )

%% -----------------------------------------------------------------------------
% APPLY WINDOWED FUNCTION TO A MOTION CAPTURE TIME SERIES

addpath(genpath('D:\Octave code\'));
wii_filename = 'square_rest_circle_period4s.wii';
fid_wii = fopen(wii_filename,'rt');
wii_headers = textscan(fid_wii,'%s',2,'Delimiter','\n');
wii_sample_frequency = wii_headers{1}{1};
wii_data = textscan(fid_wii,'%f %f %f %f %f %f %f');

wii_accel_data = wii_data(1:3);
amount_wii_dimensions = length(wii_accel_data);
wii_accel_data = cell2mat(wii_accel_data');
length_wii_timeseries = length(wii_accel_data) / amount_wii_dimensions; 
wii_accel_data = reshape(wii_accel_data,amount_wii_dimensions, length_wii_timeseries);
wii_accel_data = sqrt(sum(wii_accel_data.^2)); % magnitude

mocap_data.raw = wii_accel_data(:,2:wii_sample_frequency*24); % assign and trim;
window_size = 200; % 200
mocap_data.processed = window_function(mocap_data.raw,wii_sample_frequency,'gauss_kernel',window_size);

% plot mocap data:
% Implement multiple plot as a function (m.file):
f = multiplot_structure(mocap_data); 

% ==============================================================================
%                   Windowed functions in 2 dimensions
%% -----------------------------------------------------------------------------
% MAKE DIFFERENT 2D KERNELS

% Gaussian 2D:
gauss_2D_size = 100;          % blur = 5;    edge = 3
gauss_2D_height = 9;          % blur = 36;   edge = 9
gauss_2D_voffset = 0;         % blur = 0;    edge = -1
gauss_2D_alpha_factor = pi/2; % blur = 1.47; edge = pi
gauss_2D_alpha = gauss_2D_alpha_factor*pi/gauss_2D_size;
[kernel_gauss_x, kernel_gauss_y] = meshgrid((-(gauss_2D_size-1)/2):((gauss_2D_size-1)/2), (-(gauss_2D_size-1)/2):((gauss_2D_size-1)/2));
kernels_2D.gauss = exp( -((gauss_2D_alpha.*kernel_gauss_x).^2) / 2 - ((gauss_2D_alpha.*kernel_gauss_y).^2) / 2);
kernels_2D.gauss = kernels_2D.gauss .* gauss_2D_height + gauss_2D_voffset;
%close all; figure; imagesc(kernels_2D.gauss); figure; surf(kernels_2D.gauss)
%kernels_2D.gauss

% Sombrero Mexicano:
mexihat_2D_size = 100;       % contour = 3
mexihat_2D_height = 5;       % contour = 5
mexihat_2D_alpha_factor = 1; % contour = pi*6
mexihat_2D_alpha = mexihat_2D_alpha_factor*pi/mexihat_2D_size;
[kernel_mexihat_x, kernel_mexihat_y] = meshgrid((-(mexihat_2D_size-1)/2):((mexihat_2D_size-1)/2), (-(mexihat_2D_size-1)/2):((mexihat_2D_size-1)/2));
parabola = sqrt (mexihat_2D_alpha.*kernel_mexihat_x.^2 + mexihat_2D_alpha.*kernel_mexihat_y.^2) + eps; % eps prevents division by zero
kernels_2D.mexihat = sin (parabola) ./ parabola;
kernels_2D.mexihat = kernels_2D.mexihat .* mexihat_2D_height;
%close all; figure; imagesc(kernels_2D.mexihat); figure; surf(kernels_2D.mexihat); 
%kernels_2D.mexihat

% visualise 2D kernels:
figure;
subplot(2,2,1);
imagesc(kernels_2D.gauss);
axis square
subplot(2,2,2);
imagesc(kernels_2D.mexihat);
axis square
subplot(2,2,3);
surf(kernel_gauss_x, kernel_gauss_y, kernels_2D.gauss);
axis square
subplot(2,2,4);
surf(kernel_mexihat_x, kernel_mexihat_y, kernels_2D.mexihat);
axis square

%% -----------------------------------------------------------------------------
% APPLY FILTER TO AN IMAGE

image_filename = 'D:\Octave code\juan_2005_mexico_GREY.jpg';
images.raw = imread(image_filename);

% View original image:
figure;
imshow(images.raw)
title('raw')

% convolve original image with kernels using conv2:
images.gauss = conv2(images.raw,kernels_2D.gauss,'same');
images.gauss = images.gauss - min(min(images.gauss)); % shift
images.gauss = images.gauss / max(max(images.gauss)); % normalise it, don't criticise it

images.mexihat = conv2(images.raw,kernels_2D.mexihat,'same');
images.mexihat = images.mexihat - min(min(images.mexihat)); % shift
images.mexihat = images.mexihat / max(max(images.mexihat)); % normalise it, don't criticise it

images.gaussmexihat = conv2(images.gauss,kernels_2D.mexihat,'same');
images.gaussmexihat = images.gaussmexihat - min(min(images.gaussmexihat)); % shift
images.gaussmexihat = images.gaussmexihat / max(max(images.gaussmexihat)); % normalise it, don't criticise it

images.mexigauss = conv2(images.mexihat,kernels_2D.gauss,'same');
images.mexigauss = images.mexigauss - min(min(images.mexigauss)); % shift
images.mexigauss = images.mexigauss / max(max(images.mexigauss)); % normalise it, don't criticise it

% View original and filtered:
figure;
imshow(images.raw)
title('raw')
figure;
imshow(images.gauss);
title('gauss')
figure;
imshow(images.mexihat);
title('mexihat')
figure;
imshow(images.gaussmexihat);
title('gaussmexihat')
imshow(images.mexigauss);
title('mexigauss')

%% -----------------------------------------------------------------------------
% GET NOVELTY SCORE OF A 1D-TIME SERIES

% make a similarity matrix of 1D-TIME SERIES (artificial time series, mocap and audio):
nov_tseries.raw = mocap_data.raw;
nov_tseries.proc = mocap_data.processed;
nov_tseries.proc = nov_tseries.proc - min(min(nov_tseries.proc));
nov_tseries.proc = nov_tseries.proc / max(max(nov_tseries.proc));
nov_simmat = squareform(pdist(nov_tseries.proc'));
close all; figure; imagesc(nov_simmat)
title('similarity matrix')

% make a Gaussian Checkerboard Kernel:
gausscb_2D_size = 100; % 100
gausscb_2D_height = 1; % 1
gausscb_2D_alpha_factor = 2; % 2
gausscb_2D_alpha = gausscb_2D_alpha_factor*pi/gausscb_2D_size;
[kernel_gausscb_x, kernel_gausscb_y] = meshgrid((-(gausscb_2D_size-1)/2):((gausscb_2D_size-1)/2), (-(gausscb_2D_size-1)/2):((gausscb_2D_size-1)/2));
kernels_2D.gausscb = exp( -((gausscb_2D_alpha.*kernel_gausscb_x).^2) / 2 - ((gausscb_2D_alpha.*kernel_gausscb_y).^2) / 2);
kernels_2D.gausscb = kernels_2D.gausscb .* gausscb_2D_height;
checkerboard = kron([-1,1;1,-1],ones(gausscb_2D_size/2));
figure; imagesc(checkerboard);
title('checkerboard')
kernels_2D.gausscb = kernels_2D.gausscb .* checkerboard;
figure; imagesc(kernels_2D.gausscb)
title('Gaussian checkerboard')


% add inner-margin-average padding to beginning and ending of similarity matrix's diagonal
length_nov_simmat_avpad = length(nov_tseries.proc) + gausscb_2D_size;
margin_nov_simmat_avpad = fix(gausscb_2D_size/2);
nov_simmat_avpad = zeros(length_nov_simmat_avpad,length_nov_simmat_avpad);
beginning_average = nov_simmat(1:margin_nov_simmat_avpad,1:margin_nov_simmat_avpad);
beginning_average = mean(beginning_average(:));
nov_simmat_avpad(1:length_nov_simmat_avpad,1:length_nov_simmat_avpad) = beginning_average;
ending_average = nov_simmat( length(nov_tseries.proc)-margin_nov_simmat_avpad:length(nov_tseries.proc), length(nov_tseries.proc)-margin_nov_simmat_avpad:length(nov_tseries.proc) );
ending_average = mean(ending_average(:));
nov_simmat_avpad(1:length_nov_simmat_avpad,1:length_nov_simmat_avpad) = ending_average;
nov_simmat_avpad( margin_nov_simmat_avpad+1:length(nov_tseries.proc)+margin_nov_simmat_avpad, margin_nov_simmat_avpad+1:length(nov_tseries.proc)+margin_nov_simmat_avpad ) = nov_simmat;
figure; imagesc(nov_simmat_avpad)
title('padded similarity matrix')

% convolve Gaussian Checkerboard kernel along diagonal of similarity matrix:
hop = 1;
nov_tseries.novelty_padded = zeros(1,length_nov_simmat_avpad); % init
for i = [1:hop:length(nov_tseries.proc)]
  window_start = i;
  window_end = window_start + gausscb_2D_size - 1;
  function_window = zeros(gausscb_2D_size,gausscb_2D_size);
  function_window = nov_simmat_avpad(window_start:window_end,window_start:window_end);
  function_output = sum(sum(function_window.*kernels_2D.gausscb));
  nov_tseries.novelty_padded(1,i + margin_nov_simmat_avpad) = function_output;
end

% output of same length (get rid of outer margins):
nov_tseries.novelty_same = nov_tseries.novelty_padded(margin_nov_simmat_avpad +1 : end - margin_nov_simmat_avpad ); 
nov_tseries.novelty_same = nov_tseries.novelty_same - min(min(nov_tseries.novelty_same));
nov_tseries.novelty_same = nov_tseries.novelty_same / max(max(nov_tseries.novelty_same));

% select novelty peaks over a threshold:
peak_threshold = 0.2;
peaks_boolind = [0, (diff(sign(diff(nov_tseries.novelty_same))) == -2), 0];
nov_peaks = nov_tseries.novelty_same .* peaks_boolind;
nov_peaks = (nov_peaks >= peak_threshold);
nov_peaks = nov_tseries.novelty_same .* nov_peaks;
nov_peaks(nov_peaks==0) = NaN;

% visualise stages of computing novelty peaks:
figure;
subplot(3,1,1) 
plot(nov_tseries.raw)
title('novelty raw')
subplot(3,1,2) 
plot(nov_tseries.proc)
title('novelty processed')
subplot(3,1,3) 
plot(nov_tseries.novelty_same)
hold on
plot(nov_peaks,'.r','markersize',10)
title('novelty same and peaks')

peaks_indexes = find( nov_peaks > 0 )
