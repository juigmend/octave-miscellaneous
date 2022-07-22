%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                  LINEAR FIT                                  %
%                                                                              %
%                                                        Juan Ignacio Mendoza  %
%                                                            doctoral student  %
%                                   Music Department, University of Jyväskylä  %
%                                                               January, 2016  %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested in Octave 4.

% ==============================================================================
% Initialisation:
clc
clear all
close all
% ------------------------------------------------------------------------------
% Description:

% Two methods are shown to model a labelled dataset.
% In fancy words these are known as techniques of "machine learning".

% ------------------------------------------------------------------------------
% Example:

x = [1 1 1 2 2 2 3 3 3 4 4 4]'; % labels of data points (Independent variable)
y = [1 2 3 3 6 8 4 7 9 5 8 9]'; % data points (Dependent variable)

% Polynomial fit of averages will not work unless:
% x should be consecutive natural numbers starting from 1
% amount of elements in x should be the same for all labels

% ------------------------------------------------------------------------------
% Set parameters:

poly_deg = 4; % polynomial fit degree (1 is a simple line, more is a curve)

% ------------------------------------------------------------------------------
% plot data:

scrsz = get(0,'ScreenSize'); % get screen size
fig = figure('Position',[scrsz(3)/2.2, scrsz(4)/12, scrsz(3)/2 ,scrsz(4)/2]); % set figure position and size
plot(x,y,'.','MarkerSize',14,'Color', [0, 0.8, 0.2]); % green dots
hold on

% ------------------------------------------------------------------------------
% Plot polynomial fit of averages:

% make matrix of data points:

data_matrix = zeros(size(y,1)/max(x),max(x));

for class_counter = 1:max(x)

row_counter = 1;

   for y_counter = 1:size(y,1)
      
      if x(y_counter) == class_counter
         
         data_matrix(row_counter,class_counter) = y(y_counter);
         row_counter = row_counter + 1;
         
      end 
   end
end

means = mean(data_matrix,1); % averages
indep_var_values = [1:max(x)];

p = polyfit(indep_var_values,means,poly_deg); % polynomial fit coefficients
regline = polyval(p,indep_var_values); % evaluation of coefficients

plot(indep_var_values,regline,'linestyle','-','linewidth',10,...
      'Color', [1, 0.7, 0.1]) % plot polynomial fit (orange contiuous thick line)

% ------------------------------------------------------------------------------
% Plot linear regression using normal equation:

% make x and y:

counter_vector = 1;

x = (zeros(1, ( size(data_matrix,1) * size(data_matrix,2) ) ))';
y = x;

for index_x = 1:size(data_matrix,2)
   
   for index_y = 1:size(data_matrix,1)
      
      x(counter_vector,1) = index_x;
      y(counter_vector,1) = data_matrix(index_y,index_x)';
      counter_vector = counter_vector + 1;
      
   end
end

X = [ones(size(x,1), 1), x];
theta = ( pinv(X'*X) ) * X' * y;

plot( (X(:,2)), (X * theta), 'linestyle','--','linewidth',4,...
      'Color', [0.2, 0.2, 1] ) % plot linear regression (blue dashed thin line)

% ------------------------------------------------------------------------------
% Format plot:

axis([0.5 4.5 0 10])
set(gca,'xtick',[1:max(x)])
xlabel('x','fontsize',18)
ylabel('y','fontsize',18)

legend('Training data',['polyfit deg. = ',num2str(poly_deg)],'normal fit','location','SouthEast') 
  
% ------------------------------------------------------------------------------
% Save the figure:

%saveas(fig,strcat('SIMPLE_LINEAR_FIT_demo','.png'))
%saveas(fig,strcat('SIMPLE_LINEAR_FIT_demo','.pdf'))