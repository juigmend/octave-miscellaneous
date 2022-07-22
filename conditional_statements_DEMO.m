%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                            OCTAVE / MATLAB TUTORIAL                          %
%                             CONDITIONAL STATEMENTS                           %
%                                                                              %
%                                 October, 2017                                %
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

% A conditional statement checks if a condition is met, then do (or not do) 
% something else.

%% -----------------------------------------------------------------------------
% Example 1:

a = 0;  % a variable is initialised with a value
b = 10; % another variable is initialised with a value

% Check if a and b are equal, if the result is 1 ('true' or 'yes') then do 
% everything until 'end':
if a == b 

  display('*** a and b are equal ***')

end


% Check if a and b are not equal, if the result is 1 ('true' or 'yes') then do 
% everything until 'end':
if a ~= b
 
  display('*** a and b are NOT equal ***')

end


%% -----------------------------------------------------------------------------
% Example 2:

% The previous example can be written shorter as follows:

if a == b % Check if a and b are equal, if they are then do what follows

  display('*** a and b are equal ***')

else % In any other case ('else'), do what follows
 
  display('*** a and b are NOT equal ***')

end


%% -----------------------------------------------------------------------------
% Example 3:

% A more detailed comparison between a and b:


if a == b % Check if a and b are equal, if they are then do what follows

  display('*** a and b are equal ***')

elseif a > b  % If a is greater than b, do what follows
 
  display('*** a is greater than b ***')
  
elseif a < b  % If a is less than b, do what follows
 
  display('*** a is less than b ***')

end


%% -----------------------------------------------------------------------------
% Example 4:

% A more interesting example, in which we search if a name exists in a database.
% The function 'strcmp' is used to check if two character strings (i.e., text)
% are equal.

% Make database:
names = {'Sini','Tobias','Sumeeta', 'Alvaro', 'Anna', 'Deniz', 'Lauri', 'Kirsti', 'Riku', 'Anastasios', 'Manu', 'Kendra', 'Joonas', 'Frankie', 'Farshad', 'Saana', 'Santeri', 'Pasi', 'Juan', 'Juan'};

query = 'Sumeeta'; % <--- name to search

% 4a)
% Go through the database and compare each entry with the query:
for i = 1:length(names)
  
  if strcmp(names{i},query)
    
    display('The name exists')
    
  end

end

%% 4b)
% The following code will also tell the location of the name in the database,
% the locations of the name if it appears more than once and if the name doesn't
% exist:

query = 'Deniz';  % <--- name to search
% query = 'Juan';   % <--- name to search
% query = 'Donald'; % <--- name to search

boolean_index = strcmp(names,query) % logical indexing instead of a loop 

if sum(boolean_index) == 0 % notice the use of the 'sum' function
  
  display('The name does not exist in the database :(')

elseif sum(boolean_index) == 1
  
  location = find(boolean_index);
  display(['The name exists at location ',num2str(location)])

elseif sum(boolean_index) > 1 % this is for a name that repeats

  location = find(boolean_index);
  display(['The name exists at locations ',num2str(location)])

end