%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                            OCTAVE / MATLAB TUTORIAL                          %
%                            FIND VALUES IN AN ARRAY                           %
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

% Below are some methods to find values in an array.
% Run each line separately.

%% -----------------------------------------------------------------------------
% MATRICES AND CELLS

% First, we need to be aware of two types of arrays: matrix/vector and cell.

% Example of a Vector:

v = [7 8 9 10 22 134]
v(1)     % this returns the first element in the vector (7)
v(2:3)   % elements 2 to 3 (8, 9)
v(3:end) % elements 3 till the end (9 10 22 134)

% Example of a Matrix:

m = [11 12 13; 14 15 16; 17 18 19]
m(2,3) % this returns the element of m at row 2 and column 3

% Vector-indexing of a matrix:
m(:)   % all elements of a matrix as a vector, as if columns were concatenated
m(3:7) % elements 3 to 7, of the matrix as a vector

% A variable can be of the type 'char' (i.e., character):
my_variable  = 'a string of characters'

% However, see what happens when several words are stored as if each was a 
% string and each string an element of a vector:
failed_character_vector_1 = ['a','collection','of','strings']

% Could we store them separately by indexing explictly each word?
failed_mixed_vector_2(1) = 'a';
failed_mixed_vector_2(2) = 'collection';
failed_mixed_vector_2(3) = 'of';
failed_mixed_vector_2(4) = 'strings'

% A cell is a type of array that can hold several values of different dimensions, 
% for example several matrices or several separate strings of alphanumeric 
% characters (i.e., text). 
% Cells can be nested (one or more cells in a containing  cell).
% Cells are written between curly brackets.

% Examples of Cells:

c1 = {[1 2 3 4],[50; 60; 70; 80]} % a cell containing a row vector and a column vector

c2 = {'cat','dog'} % a cell containing character strings

c3 = {[74 56 128], [10,20; 30, 40], 'cat'} % a cell containing a mixture

% Retrieve values from a cell array:

c3(1) % returns a cell containing vector [74 56 128]
c3(2) % returns a cell containing matrix [10,20; 30, 40]
c3(3) % returns a cell containing string 'cat'

c3{1} % returns vector [74 56 128]
c3{2} % returns matrix [10,20; 30, 40]
c3{3} % returns string 'cat'

c3{1}(3)   % returns 128
c3{2}(2,1) % returns 30
c3{3}(2,1) % returns an error. Why? Think about it.

% Hence, this should work nicely:
successful_mixed_cell{1} = 'a';
successful_mixed_cell{2} = 'collection';
successful_mixed_cell{3} = 'of';
successful_mixed_cell{4} = 'strings'

%% -----------------------------------------------------------------------------
% FIND MAXIMUM AND MINIMUM OF AN ARRAY

% The functions 'max' and 'min' find maximum or minimum in a matrix or vector.
% In a matrix, they will return values column-wise, unless stated otherwise
% (see 'help max' and 'help min').

max(v)      % returns 9
min(m)      % returns [11 12 13]
min(min(m)) % returns 11
min(m(:))   % returns 11

%% -----------------------------------------------------------------------------
% FIND BY LOGICAL INDEXING

% Logical data is comprised by 'Boolean' values 'true' (or 'yes') 
% and 'false' (or 'no'). They are usually represented as 1 and 0, respectively. 
% George Boole (1815 - 1864) was an English mathematician.

a = [1, 2, 5, 8 , 99, -4, 6, 22, - 200.9]

b = 8;

% c is a logical vector of same dimensions as a,
% with a 1 at the location of values in b
c = a == b 

% Other examples:

a < 0   % all negative values
a >= 22 % all values greater or equal than 22

% The index of a boolean value can be found with the function 'find':

lonely_one = [0 0 0 0 1 0 0]

index_of_lonely_one = find(lonely_one)

% Thus, the following will return the index of b in a:
find(a == b)

%% -----------------------------------------------------------------------------

% In what follows, an example in which an input number is searched in a database
% where each number is linked to a tag (i.e., a text). If the input number 
% matches a number in the database, the program will return its associated tag.
% Otherwise, the program will return the tag associated to the closest number.
% To make it more interesting, a 'real life problem' is proposed.

% Problem:
% At the finish of the marathon of Jyv??skyl?? each runner's name was recorded 
% in a database, along with their time of completion (in seconds).
% We want a computer program that accepts as input a time and returns the name
% of the runner closest to that time and the recorded time for that name.

runners = {'Sini','Tobias','Sumeeta', 'Alvaro', 'Anna', 'Deniz', 'Lauri', 'Kirsti', 'Riku', 'Anastasios', 'Manu', 'Kendra', 'Joonas', 'Frankie', 'Farshad', 'Saana', 'Santeri', 'Pasi', 'Juan'} 
  times = [ 1704 , 1746   , 1802    ,  1850   ,  1853 ,  1855  ,  1930  ,  1932   ,  2010 ,  2012       ,  2039 ,  2045   ,  2048   ,  2057    ,  2058    ,  2059  ,  3207    ,  3209 ,  3210 ]      

  query = 3208; % <--- input time

difference = abs(query - times) % absolute difference between query and all recorded times
closest_value = min(difference) % minimum value of the differences (the closest to the query)
closest_boolean = difference == closest_value
closest_index = find(closest_boolean)
associated_name = runners{closest_index}
associated_time = times(closest_index)

display(['Runner with a time closest to the query: ',associated_name])
display(['time: ',num2str(associated_time(1)),' (s)'])

% Note that if the query is exactly between two values in the database, the
% returned values will be only one time and name.
% Can you figure out how to display both times and names?
