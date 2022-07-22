%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                            OCTAVE / MATLAB TUTORIAL                          %
%                                     LOOPS                                    %
%                                                                              %
%                                September, 2017                               %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INSTRUCTIONS:

% This tutorial is written in code for GNU Octave or Matlab.
% Read the comments and run each cell of code (commands) separately.
% The comments are preceded by a percentage sign (%).
% A cell is delimited by two percentage signs (%%). 
% In Matlab A cell can be run by hitting the keys "command" and "enter" at the 
% same time.
%%

% This is a cell.

%% =============================================================================

% A loop is  a process repeating several times, where each time some of its 
% variables change.

%...............................................................................

% First, familiarity with basic logical operators is needed.
% Logical operators compare two values.
% They can return two values:
% 1 is 'true' (or 'yes')
% 0 is 'false' (or 'no')
% They can be thought of as a question.

a == b  % are a and b equal?
a ~= b  % are a and b different? (Matlab)
a != b  % are a and b different? (Octave)
a > b   % is a greater than b?
a < b   % is a less than b?
a >= b  % is a greater or equal than b?
a <= b  % is a less or equal than b?

%% -----------------------------------------------------------------------------
% LOOP METHOD 1: 'while' loop
clc

a = 0;  % a variable is initialised with a value
b = 10; % another variable is initialised with a value

% While a is less or equal than b do everything below intil 'end':
while a <= b 
    
    a = a + 1 % this will add 1 to a
        
end

% This loop simply counts from 1 to 11
% Note that when 'while' checks that a is equal to b, it will perform the process
% for one last time.

%% .............................................................................
% A more interesting example of a 'while' loop:
clc

lonely_one = [0 0 0 0 0 0 0 0 1 0 0 0 0 ]; % a one in the middle of zeroes

% find the lonely one:
index = 1;
a = 0;
while a == 0 
    
    index = index + 1;
    a = lonely_one(index) == 1;
        
end

display(['the index of the lonely one is ',num2str(index)])

% Try changing the location of the one in the 'lonely_one' vector.

%% -----------------------------------------------------------------------------
% LOOP METHOD 2: 'for' loop
% comment/uncomment one assignment to variable d, try them all:
clc

d = [1, 2, 3]; 
% d = [1, 3, 5];
% d = 1:10;

% Do the process below until 'end' only for the values in variable d,
% which will be assigned one by one to c:
for c = d 

    c % just print c

end

%% ............................................................................. 
% A more interesting example of a 'for' loop:
clc

%%

e = [ 1 2 1 3 1 2 1 0 0 0 0 0 0 0 ]; % a short signal
plot(e)

%%

e_1 = repmat(e,1,10); % a longer signal is made by repeating the short signal
plot(e_1)

%%
% The loop below will do this:
%   Compute the average of the first 4 elements of e_1 (this is the elements at indexes 1 to 4).
%   Then compute the average of elements of e_1 at indexes 2 to 5. 
%   Then compute the average of elements of e_1 at indexes 3 to 6. 
%   Continue doing this process until reaching the end of e_1 (indexes end-3 to end).
%   Store all the results as consecutive values in a vector.

window_length = 4; % every time we compute an average we take a 'window' of 4 consecutive elements

% note that the 'counting' variable can have any name (this one is a tribute to a famous musician)
for snoopdogg = 1:( length(e_1) - window_length ) % this is the amount of repetitions of the process (think about it)
    
    window_start = snoopdogg;
    window_end = window_start + window_length;
    result(snoopdogg) = mean( e_1(window_start:window_end) ); % here the average (A.K.A. 'mean') is computed upon the window

end

plot(result)

% The 'for' loop above is the classical "moving average filter".
% It is widely used to smooth a signal, as you can see in the last plot.
% Try using different values for 'window_length'.
