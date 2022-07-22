function output = window_function(varargin);

%  output = window_function(time_series,[empty],window_function, window_size, hop);
%
%  input:
%  time_series (matrix where columns are channels and rows are samples)
%  window_function = {'mean', 'mean_power', 'square_kernel', 'triangle_kernel', 'gauss_kernel', 'fft'} 
%  window_size (samples)
%  hop (samples, optional, default = 1)
%  extension = {'same', 'average_pad', 'valid'} (samples, optional, default = 'same')
%
%  note:
%  time_series rows are dimensions and columns are samples (observations)
%  padding is always the average of the inner margin
%
%  July,2017
%  Juan Ignacio Mendoza
%  University of Jyväskylä

kernel_switch = 0;
kernel = 0;

if length(varargin) == 4
  varargin{5} = 1;
end

if length(varargin) == 5
  varargin{6} = 'same';
end

time_series = varargin{1};
window_size = varargin{4};
hop = varargin{5}; 
time_series_size = size(time_series);

if xor(strcmp(varargin{3},'square_kernel'), strcmp(varargin{3},'triangle_kernel'))
  kernel = ones(1,window_size);
  kernel_switch = 1;
end

if strcmp(varargin{3}, 'triangle_kernel')
  kernel(1:fix(window_size/2)) = 0 : 1 / ( fix(window_size/2) - 1 ) : 1 ;
  kernel(fix(window_size/2)+1:end) = flip(kernel(1:fix(window_size/2))) ;
  kernel_switch = 1;
end

if strcmp(varargin{3}, 'gauss_kernel')
  gauss_alpha = 2*pi/window_size;
  gauss_height = (sqrt(2*pi))/window_size;
  x = linspace( -(window_size-1)/2, (window_size-1)/2, window_size );
  kernel = gauss_height * ( exp( -((gauss_alpha *x).^2) / 2 ) );
  kernel_switch = 1;
end

% add zero-pad margins:
extended_time_series_zeropad = zeros(time_series_size(1),time_series_size(2) + window_size);
window_margin = fix(window_size/2); 
for i = 1:time_series_size(1)
  extended_time_series_zeropad( i,window_margin + 1 : (time_series_size(2) + window_margin)) = time_series(i,:);
end
extended_time_series_padded = extended_time_series_zeropad; 

% add inner-margin-average padding:
extended_time_series_meanpad = extended_time_series_zeropad; % init the mean pad :(
for i = 1:time_series_size(1)
  extended_time_series_meanpad(i,1:window_margin) = mean(time_series(i,1:window_margin));
  extended_time_series_meanpad(i,window_margin+time_series_size(2):end) = mean(time_series(i,time_series_size(2)-window_margin+1:end));
end
extended_time_series_padded = extended_time_series_meanpad; 
amt_hopwindows = fix(time_series_size(2)/hop);

wbar = waitbar(0,'Busy.');
wbarfrac = 1 / (time_series_size(1) * (amt_hopwindows+1));
wbarcurr = 0;

% %----------------------------------------------------------------------------
% % NON-OPTIMAL SOLUTION:
%
%hopslide_padded = zeros(size(extended_time_series_padded)); % init
%
%for i_1 = 1:time_series_size(1) % dimensions
%  
%  for i = [1:hop:time_series_size(2)] % sliding window
%    
%    window_start = i;
%    window_end = window_start + window_size - 1;
%    function_window = zeros(1,window_size);
%    function_window = extended_time_series_padded(i_1,window_start:window_end);
%    
%    if strcmp(varargin{3},'mean')
%      function_output = mean(function_window);
%    end
%    
%    if strcmp(varargin{3},'mean_power')
%      function_output = (mean(function_window'))^2;
%    end
%    
%    if kernel_switch == 1;
%      function_output = sum(function_window.*kernel);
%    end
%    
%    if strcmp(varargin{3},'fft')
%      function_output = fft(function_window);
%    end
%    
%    hopslide_padded(i_1,i + window_margin) = function_output;
%    
%    wbarcurr = wbarcurr + wbarfrac;
%    waitbar(wbarcurr);
%  end
%end
%
% %----------------------------------------------------------------------------
% % GOOD SOLUTION:
%
%hopslide_padded = zeros(size(extended_time_series_padded)); % init
%
%if strcmp(varargin{3},'mean')
%  function_handle = @mean;
%end
%
%if strcmp(varargin{3},'mean_power')
%  function_handle = @(function_window) (mean(function_window'))^2;
%end
%
%if kernel_switch == 1;
%  function_handle = @(function_window) sum(function_window.*kernel);
%end
%
%if strcmp(varargin{3},'fft')
%  function_handle = @fft;
%end
%
%for i_1 = 1:time_series_size(1) % dimensions
%  
%  for i = [1:hop:time_series_size(2)] % sliding window
%    
%    window_start = i;
%    window_end = window_start + window_size - 1;
%    function_window = zeros(1,window_size);
%    function_window = extended_time_series_padded(i_1,window_start:window_end);
%    
%    function_output = function_handle(function_window); % General Window Function   
%    
%    hopslide_padded(i_1,i + window_margin) = function_output;
%    
%    wbarcurr = wbarcurr + wbarfrac;
%    waitbar(wbarcurr);
%  end
%end

% %----------------------------------------------------------------------------
% % BETTER SOLUTION:

if strcmp(varargin{3},'mean')
  function_handle = @mean;
end

if strcmp(varargin{3},'mean_power')
  function_handle = @(function_window) (mean(function_window'))^2;
end

if kernel_switch == 1;
  function_handle = @(function_window) sum(function_window.*kernel);
end

if strcmp(varargin{3},'fft')
  function_handle = @fft;
end

hopwindows_margin = fix(window_margin/hop); 
hopslide_padded = zeros( 1 , amt_hopwindows + window_size ); % init

for i_1 = 1:time_series_size(1) % dimensions
  
  for i = [1:amt_hopwindows] % sliding window
    window_start = i*hop - hop + 1;
    window_end = window_start + window_size - 1;  
    hopslide_padded(i_1,i + hopwindows_margin) = function_handle(extended_time_series_padded(i_1,window_start:window_end)); % General Window Function
    wbarcurr = wbarcurr + wbarfrac;
    waitbar(wbarcurr);
  end
end

if strcmp(varargin{6},'average_pad')
  output = hopslide_padded;
end
 
if xor(strcmp(varargin{6},'same'), strcmp(varargin{6},'valid'))
  output = hopslide_padded(:,hopwindows_margin +1 : end - hopwindows_margin);
end

if strcmp(varargin{6},'valid')
  output = output( hopwindows_margin +1 : end - hopwindows_margin ); 
end

close(wbar)

end