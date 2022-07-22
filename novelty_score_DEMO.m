%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                 NOVELTY SCORE                                %
%                                                                              %
%                                  July, 2018                                  %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyväskylä                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tested on Octave 4

%===============================================================================

% DESCRIPTION:

%   Produce a signal with some randomness, then a self-similarity matrix of it.
%   The matrix is bounded, to reduce computing time.
%   Then convolve a Gaussian-tapered checkerboard kernel along the diagonal of 
%   the self-similarity matrix. The result is a novelty score.
%   (c.f. Foote & Cooper, 2003)

% INSTRUCTIONS:

%   1) Modify the values indicated by an arrow (<---)
%   2) Run the script and enjoy.

%===============================================================================

clc
clear
close all

%===============================================================================

   signal_type ='sqr'; % <--- signal carrier type ('sin' or 'sqr')
 signal_length = 99;   % <--- signal length
 signal_period = 33;   % <--- signal period
     bandwidth = 11;   % <--- gaussian-tapered checkerboard kernel bandwidth
       full_sw = 0;    % <--- 1 = upper and lower triangles; 0 = upper triangle

%-------------------------------------------------------------------------------

if strcmp(signal_type,'sin')
  x = [1:signal_length]; % grid
  signal = sin(x*pi*2/(signal_period) + randn(1)/8) + (sin(x*pi/4*randn(1)*2)/6+randn(1)/100);
elseif strcmp(signal_type,'sqr')
  signal = ones(1,signal_period);
  signal(ceil(signal_period/2):end) = -1;
  signal = repmat(signal,1,ceil(signal_length/signal_period));
  signal = signal(1:signal_length) + randn(1,signal_length)/4;
end
signal = signal .* [1:signal_length]/signal_length/2; 

selfsim_matrix = zeros(signal_length); % initialise self-similarity matrix

% boundary's ending indexes:
boundary = bandwidth - 1;
end_bounds = [1:signal_length]; % increasing indexes
end_bounds(end-boundary+1:end) = ...
  end_bounds(end-boundary+1:end) - [1:boundary]; % at the end the available space shrinks
end_bounds = end_bounds + boundary; % shift for length of boundary

% compute distances;
for row = 1:signal_length
  for col = row+1:end_bounds(row) % start at diagonal + 1; end at corresponding index
    selfsim_matrix(row,col) = abs( signal(row) - signal(col) ); % Euclidean distance
  end
end

if full_sw
  selfsim_matrix = selfsim_matrix + flipud(rot90(selfsim_matrix)); % upper and lower triangles
end

% Gaussian Checkerboard Kernel:
gauss_2D_size = round(signal_length/2);
if rem(gauss_2D_size,2)
  gauss_2D_size = gauss_2D_size + 1; % has to be even to make the checkerboard
end
gauss_2D_height = 1/(sqrt(2*pi)*bandwidth);
gauss_2D_variance = bandwidth^2;
[kernel_gauss_x, kernel_gauss_y] = meshgrid((-(gauss_2D_size-1)/2):((gauss_2D_size-1)/2),...
  (-(gauss_2D_size-1)/2):((gauss_2D_size-1)/2));
kernels_2D.gauss = exp( -(kernel_gauss_x.^2) / (2*gauss_2D_variance) ...
  - (kernel_gauss_y.^2) / (2*gauss_2D_variance));
kernels_2D.gauss = kernels_2D.gauss .* gauss_2D_height;
kernels_2D.gausscb = kernels_2D.gauss .* kron([-1,1;1,-1],ones(gauss_2D_size/2));

% add inner-margin-average padding to beginning and ending of similarity matrix's diagonal
length_selfsim_matrix_avpad = length(signal) + gauss_2D_size;
margin_selfsim_matrix_avpad = fix(gauss_2D_size/2);
selfsim_matrix_avpad = zeros(length_selfsim_matrix_avpad,length_selfsim_matrix_avpad);
beginning_average = selfsim_matrix(1:margin_selfsim_matrix_avpad,...
  1:margin_selfsim_matrix_avpad);
beginning_average = mean(beginning_average(:));
selfsim_matrix_avpad(1:length_selfsim_matrix_avpad,1:length_selfsim_matrix_avpad) = ...
  beginning_average;
ending_average = selfsim_matrix( length(signal) - margin_selfsim_matrix_avpad : length(signal),...
  length(signal) - margin_selfsim_matrix_avpad : length(signal) );
ending_average = mean(ending_average(:));
selfsim_matrix_avpad(1:length_selfsim_matrix_avpad,1:length_selfsim_matrix_avpad) = ...
  ending_average;
selfsim_matrix_avpad( margin_selfsim_matrix_avpad+1:length(signal)+margin_selfsim_matrix_avpad,...
  margin_selfsim_matrix_avpad+1:length(signal)+margin_selfsim_matrix_avpad ) = ...
  selfsim_matrix;

tic
% convolve Gaussian Checkerboard kernel along diagonal of similarity matrix:
nov_tseries.novelty_padded = zeros(1,length_selfsim_matrix_avpad); % init
for i = [1:length(signal)]
  window_start = i;
  window_end = window_start + gauss_2D_size - 1;
  function_window = zeros(gauss_2D_size,gauss_2D_size);
  function_window = selfsim_matrix_avpad(window_start:window_end,window_start:window_end);
  function_output = sum(sum(function_window.*kernels_2D.gausscb));
  nov_tseries.novelty_padded(1,i + margin_selfsim_matrix_avpad) = function_output;
end
toc

% output of same length (get rid of outer margins):
nov_tseries.novelty_same = nov_tseries.novelty_padded(margin_selfsim_matrix_avpad + 1 :...
  end - margin_selfsim_matrix_avpad );
 
% normalise it, don't criticise it: 
nov_tseries.novelty_same = nov_tseries.novelty_same - min(min(nov_tseries.novelty_same));
nov_tseries.novelty_same = nov_tseries.novelty_same / max(max(nov_tseries.novelty_same));

% visualise:
figure;
subplot(2,2,1)
imagesc(kernels_2D.gauss); 
title('gaussian bell heatmap')
subplot(2,2,2)
surf(kernels_2D.gauss)
title('gaussian bell surface')
subplot(2,2,3)
imagesc(kernels_2D.gausscb)
title('gaussian checkerboard heatmap')
subplot(2,2,4)
surf(kernels_2D.gausscb)
title('gaussian checkerboard surface')
figure;
subplot(6,1,1)
plot(signal)
title('signal')
subplot(6,1,2:5)
imagesc(selfsim_matrix)
title('self-similarity matrix')
subplot(6,1,6)
plot(nov_tseries.novelty_same)
title('novelty score')