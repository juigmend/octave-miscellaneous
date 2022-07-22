%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                              Gaussian Function                               %
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

% The Gaussian Function is typically used as a measure of normal distribution in
% probability. 

% It is also used as a smoothing function, by using it as a window that is
% convolved with an input signal.

% The following code allows to explore the workings of the Gaussian Function in
% one dimension. The vector w can be used as a kernel for smoothing.

% Mind that the Gaussian Function is continuous and the method presented here is
% discrete. Therefore the computations have a small error.

% ------------------------------------------------------------------------------
% 1 - Input Parameters
% Comment/uncomment options 1.a or 1.b

%% 1.a - Height = 1 for any bandwidth:
%bandwidth = 15
%alpha = 2*pi/bandwidth
%height = 1

% 1.b - Area = 1 for any bandwidth:
bandwidth = 15
alpha = 2*pi/bandwidth
height = (sqrt(2*pi))/bandwidth

% ------------------------------------------------------------------------------
% 2 - Compute Gaussian Function
% Comment/uncomment options 2.a or 2.b

% 2.a - use Octave's function (from the Signal package):
%pkg load all
%w = height * gaussian(bandwidth,alpha);
%plot(w,'b')

%% 2.b - embedded code:
x = linspace( -(bandwidth-1)/2, (bandwidth-1)/2, bandwidth );
w = height * ( exp( -((alpha*x).^2) / 2 ) ); 
plot(x,w,'b')

approximate_area = sum(w) % rough numerical integration
