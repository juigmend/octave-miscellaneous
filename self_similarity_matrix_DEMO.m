%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                             SELF-SIMILARITY MATRIX                           %
%                                                                              %
%                                  June, 2018                                  %
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
      boundary = 28;   % <--- amount of distances from each value 
       full_sw = 1;    % <--- 1 = upper and lower triangles; 0 = upper triangle

%-------------------------------------------------------------------------------

if strcmp(signal_type,'sin')
  x = [1:signal_length]; % grid
  signal = sin(x*pi*2/(signal_period) + randn(1)/8) + (sin(x*pi/4*randn(1)*2)/6+randn(1)/100);
elseif strcmp(signal_type,'sqr')
  signal = ones(1,signal_period);
  signal(ceil(signal_period/2):end) = -1;
  signal = repmat(signal,1,floor(signal_length/signal_period));
  signal = signal + randn(1,signal_length)/4;
end
signal = signal .* [1:signal_length]/signal_length/2; 

selfsim_matrix = zeros(signal_length); % initialise self-similarity matrix

% boundary's ending indexes:
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

% visualise:
subplot(5,1,1)
plot(signal)
title('signal')
subplot(5,1,2:5)
imagesc(selfsim_matrix)
title('self-similarity matrix')