function f = multiplot_structure(structure);

%  f = multiplot_structure(structure);
%
%  input:
%    structure of matrices should be organized as follows:
%      structure.(matrix_name)
%      In fields of structure rows are dimensions and columns are samples (observations)
%
%  output
%    A figure with subplots. Each subplot is a plot of a field in the 
%    structure. The title of each plot is the name of the corresponding field.
%
%  July,2017
%  Juan Ignacio Mendoza
%  University of Jyväskylä

line_width = 1;
screen_div = 2;
screen_size = get(0,'screensize');
figures_position(3:4) = screen_size(3:4) ./ screen_div ;
figures_position(1:2) = repmat(figures_position(4) ./ screen_div,1,2);

structure_names = fieldnames(structure);
amount_structure_names = length(structure_names);

f = figure('position', figures_position);

for i = 1:amount_structure_names
  subplot(amount_structure_names,1,i)
  plot(structure.(structure_names{i})','linewidth',line_width)
  xlim([1,length(structure.(structure_names{i}))])
  ylim([min(structure.(structure_names{i})), max(structure.(structure_names{i}))])
  set(gca,'ytick',([min(structure.(structure_names{i})), 0, max(structure.(structure_names{i}))]))  
  title(structure_names{i},'interpreter','none')
end


end