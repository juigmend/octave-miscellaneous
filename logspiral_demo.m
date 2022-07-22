%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                          Plot Logarithmic Spirals                            %
%                                                                              %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             Doctoral Student %
%                                   Music Department - University of Jyväskylä %
%                                                                   June, 2015 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize script:
clear
close all
clc

% ------------------------------------------------------------------------------
% INPUT VALUES:
% a and b are the spiral's arbitrary constants

a = 1; % <-------------- not very interesting, 1 is okay
b = -0.5; % <----------- controls the "soul" of the spiral, sign = direction
points = 360*1; % <----- amount of points of the line
base = 0; % <----------- base of the logarithm, use 0 for natural (Euler's)

algorithm = 3; % <------ select algorithm (procedure of computation)
% 1 = given theta calculate rho, equation 1 (blue):  rho = a*exp(b*theta)
% 2 = given rho calculate theta, equation 2 (red): theta = (1/b)*(log(rho/a)
% 3 = 1 and 2
% 4 = animation (green) where a = b = 1; rho = [1:it^2*4] and:
it = 20; % <------------ amount of pictures of the animation (careful!)
p = 0; % <-------------- pause between pictures in seconds

% Notice that equation 2 gives a coarse but full representation of the spiral 
% with few points.

% rho is radius and theta is angle

% ------------------------------------------------------------------------------
figure

% solve base:
if (base == 0) == 1
  base = exp(1);
end

% these two lines are logicals to select algorithm(s) to use:
which_algorithm = zeros(1,4); 
which_algorithm(algorithm) = 1;

% algorithm 1 (blue), given theta calculate rho: 
if (which_algorithm(1)+which_algorithm(3) ~= 1) ~= 1 % 1 OR 3
  theta = [1:points];
  theta = theta*(pi/180); % deg to rad
  rho(1:points) = a*(base.^(b*theta)); % equation 1
  rho = rho/max(rho); % scaling
  polar(theta,rho,'b-*') % plot
end

if (algorithm == 3) == 1
  hold on
end

% algorithm 2 (red), given rho calculate theta:
if (which_algorithm(2)+which_algorithm(3) ~= 1) ~= 1 % 2 OR 3
  rho2 = [1:points];
  theta2 = (1/b)*((log(rho2/a))/(log(base))); % equation 2, note the change of base
  rho2 = rho2/max(rho2); % scaling
  polar(theta2,rho2,'r-*') % plot
end

% algorithm 4 (green), animation:
if (algorithm == 4) == 1
  for i = 1:it
    rho4 = [1:it*4*i];
    theta4 = ((i^2)/i)*log(rho4); % equation 2
    rho4 = rho4/max(rho4); % scaling
    refresh
    polar(theta4,rho4,'*g'); % plot
    pause(p)
   end
end

disp 'done'





