%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                            OCTAVE / MATLAB TUTORIAL                          %
%                     MULTIPLICATION OF MATRICES AND VECTORS                   %
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
% Read the comments and run each line of code (commands) separately.
% The code is the instructions to the computer (e.g., add two numbers, store a 
% number in memory and give a name to it)
% The comments are preceded by a percentage sign (%).
% Comments are not executed by the computer; their purpose is to give complementary 
% information to a person reading the code.
% A command line usually ends with a semicolon sign (;). If the semicolon is not
% present, then the result of that line will be printed in the Command Window
% To run a command line you can select the line and copy-paste it
% to the command window, then hit the "enter" key.
% Also a line can be run by selecting it and then evaluating it (in Matlab 
% control-click the selection and click "Evaluate Selection").
% I Matlab it is also possible to isolate code as a "cell" and run cells separately.
% A cell is delimited by two percentage signs (%%). A cell can be run by hitting
% the keys "command" and "enter" at the same time.
%%

% This is a cell.

%% =============================================================================

% In this context "number" may also be called "scalar", and can be stored in memory.
% Also it can be given a nam witht he only conditin that the name does not start
% with a number. This artifact for storing numbers is called a "variable". 
% The value of a variable can be changed.

variablename = 3.14 % this is a number stored in a variable called 'variablename'
variablename = -1   % this number has replaced the previous value of the variable

% An array is a way of storing several values (e.g., numbers) as a single variable
% and is written between square brackets [].

othervariablename = [1, 2, 3]        % this is an array of 3 elements stored in a variable called 'othervariablename'
              dog = 500              % this is not an array as it only has one value
              cat = [10, 10000, 0.3] % this is an array
                a = 0

% Arithmetic operations can be performed with variables:
b = cat + dog;
c = b - 4

% A matrix is an array. Therefore it is written between square brackets.
% , (comma) separates columns and ; (semicolon) separates rows
% For ease of use the comma can be replaced by an empty space.
% Matrices can have different 'dimensions', which means the amount of
% rows and the amount of columns (rows x columns)

example_matrix_1 = [1, 2, 3; 4, 5, 6; 7, 8, 9] % 3x3 matrix with commas
example_matrix_2 = [1 2 3; 4 5 6; 7 8 9]       % 3x3 matrix without commas
example_matrix_3 = [10 20 30; 40 50 60]        % 2x3 matrix
example_matrix_4 = [10 20; 30 40; 50 60]       % 3x2 matrix 

% Example_matrix_1 and example_matrix_2 are equal. 
% To select a section of a matrix, the boundaries of the section are specified
% by their indexes of rows and columns. The indexes are written between 
% parenthesis following the name of the array.

example_matrix_1(1,2)   % the element at row 1 and column 2
example_matrix_1(1,:)   % all the elements at row 1
example_matrix_1(:,3)   % all the elements at column 3
example_matrix_1(2:3,2) % elements of rows 2 to 3 from column 2

% A vector is an array of only one column and more than one row or viceversa.

example_column_vector = [1; 2; 3]  % 3x1 vector
example_row_vector_1  = [1, 2, 3]  % 1x3 vector
example_row_vector_1  = example_matrix_1(1,:) % same as above, but extracted from example_matrix_1
example_row_vector_2  = [10 20 30] % 1x3 vector

% A transposed matrix or vector is obtained by swicthing rows for columns. 
% The operator for transposition is ' (single quote mark)

example_transposed_matrix_1 = example_matrix_1'
example_transposed_matrix_1 = [1, 4, 7; 2, 5, 8; 3, 6, 9]

% * is the multiplication operator.
% A multiplication of two scalars results in a scalar.
% ( 2*3 = 6 )
example_scalar_multiplication = 2 * 3

% Multiplication of arrays (matrices and vectors) can be done element-wise or
% by dot product.

% Element-wise multiplication takes each element of a matrix or vector and 
% multiplies it with the corresponding element of another vector, given that both
% vectors are of the same size. 
% The result is a vector of the same length of the input vectors.
% Element-wise multiplication of vectors A and B is written as A.*B

example_element_wise_matrix = example_matrix_1 .* example_transposed_matrix_1
example_element_wise_matrix = [ 1*1, 2*4, 3*7 ;
                                4*2, 5*5, 6*8 ;
                                7*3, 8*6, 9*9 ]

example_element_wise_vector = example_row_vector_1 .* example_row_vector_2
example_element_wise_vector = [ 1*10, 2*20, 3*30 ]

example_element_wise_vector = example_row_vector_1' .* example_row_vector_2'
example_element_wise_vector = [ 1*10; 2*20; 3*30 ]

% this will not work because vectors are of different sizes (dimensions):
example_element_wise_vector = example_row_vector_1' .* example_row_vector_2 

% Dot product or vector/matrix multiplication is a more complex multiplication 
% operation. The dot product of matrices or vectors A and B is written as A*B
% What follows presents different cases of dot product multiplication:

% 1) Dot product of vectors. Vectors must have the same amount of elements.
%    The first input vector should be a row and the second a column.
%    Definition: The result is a scalar which is the sum of the delement-wise 
%    multiplication of both vectors.    

example_dot_vectors = example_row_vector_1 * example_row_vector_2'
example_dot_vectors = sum(example_row_vector_1 .* example_row_vector_2)

% 2) Dot product of a matrix and a vector. The vector has to be of same length 
%    as the amount of rows of the matrix. 
%    Definition: The result is a vector with the same dimensions as
%    the input vector. Each element of the result vector is the dot product of
%    the input vector and each row of the input matrix.

example_dot_matrix_vector =  example_matrix_1 * example_row_vector_1'
example_dot_matrix_vector =  example_matrix_1 * example_column_vector

example_dot_matrix_vector = [ ( example_matrix_1(1,:) * example_column_vector ) ;
                              ( example_matrix_1(2,:) * example_column_vector ) ;
                              ( example_matrix_1(3,:) * example_column_vector ) ]

example_dot_matrix_vector = [ 1*1 + 2*2 + 3*3 ;
                              1*4 + 2*5 + 3*6 ;
                              1*7 + 2*8 + 3*9 ]
                          
% 3) Dot product of matrices. The second input matrix has to have te same amount
%    of rows as the amount of columns in the first input matrix. 
%    Definition: The result is a matrix  with the same amount of rows as the 
%    first input matrix and the same amount of columns as the second input matrix. 
%    Each element of the result matrix is the dot product of the intersected vectors.

example_dot_matrices_1 =  example_matrix_1 * example_matrix_3'

example_dot_matrices_2 =  example_matrix_1 * example_matrix_3 % this will not work

% Remember that:
example_matrix_1 = [ 1 2 3 ; 
                     4 5 6 ; 
                     7 8 9 ]

example_matrix_4 = [ 10 20 ; 
                     30 40 ; 
                     50 60 ]

example_dot_matrices_3 =  example_matrix_1 * example_matrix_4

example_dot_matrices_3 = [ ( example_matrix_1(1,:) * example_matrix_4(:,1) ) , ( example_matrix_1(1,:) * example_matrix_4(:,2) ) ;
                           ( example_matrix_1(2,:) * example_matrix_4(:,1) ) , ( example_matrix_1(2,:) * example_matrix_4(:,2) ) ;
                           ( example_matrix_1(3,:) * example_matrix_4(:,1) ) , ( example_matrix_1(3,:) * example_matrix_4(:,2) ) ]

example_dot_matrices_3 = [ ( 1*10 + 2*30 + 3*50 ) , ( 1*20 + 2*40 + 3*60 ) ;
                           ( 4*10 + 5*30 + 7*50 ) , ( 4*20 + 5*40 + 6*60 ) ;
                           ( 7*10 + 8*30 + 9*50 ) , ( 7*20 + 8*40 + 9*60 ) ]
