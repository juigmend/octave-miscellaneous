%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %
%                                 TRACE RINGS                                %
%                                                                            %
%                                 July, 2018                                 %
%                                                                            %
%                         Juan Ignacio Mendoza Garay                         %
%                              doctoral student                              %
%                Department of Music, Art and Culture Studies                %
%                           University of Jyväskylä                          %
%                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tested on Octave 4

%=============================================================================

% INSTRUCTIONS:

%   1) Modify the values indicated by an arrow (<---)
%   2) Run the script and enjoy.

%===============================================================================

clc
clear
close all

%===============================================================================


% input parameters:
   width =  20 ; % <--- Width of square frame (2*radius). Use carefully.
 rings_n =   5 ; % <--- 1 <= rings_n <= radius
hi_bound =   5 ; % <--- [0:1:rings_n]
lo_bound =   4 ; % <--- [0:1:rings_n]

% Gaussian 2D:
radius = width/2;
gauss_2D_size = width;         
gauss_2D_height = -radius;          
gauss_2D_voffset = radius;         
gauss_2D_alpha_factor = pi/2;
gauss_2D_alpha = gauss_2D_alpha_factor*pi/gauss_2D_size;
[kernel_gauss_x, kernel_gauss_y] = meshgrid((-(gauss_2D_size-1)/2):((gauss_2D_size-1)/2), (-(gauss_2D_size-1)/2):((gauss_2D_size-1)/2));
gauss_bell = exp( -((gauss_2D_alpha.*kernel_gauss_x).^2) / 2 - ((gauss_2D_alpha.*kernel_gauss_y).^2) / 2);
gauss_bell = gauss_bell .* gauss_2D_height + gauss_2D_voffset;

figure
imagesc(gauss_bell); 
figure 
surf(gauss_bell)

% trace ring:
ring_width = floor(radius/rings_n);
ring_hi = ring_width * hi_bound;
ring_lo = ring_width * lo_bound; 
ring = (gauss_bell <= ring_hi) .* (gauss_bell >= ring_lo);

figure
imagesc(ring)