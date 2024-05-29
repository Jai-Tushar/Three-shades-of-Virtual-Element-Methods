function plotsol(Element, Node, solution)
% PLOT_SOLUTION plots the vertex values of a virtual element function
% original - AUTHOR: Oliver Sutton, 2016
% Modified by: Jai Tushar, BITS-Pilani, 2019



max_n_vertices = max(cellfun(@length, Element));
padding_function = @(vertex_list) [vertex_list...
			NaN(1,max_n_vertices-length(vertex_list))];
elements = cellfun(padding_function, Element, 'UniformOutput', false);
elements = vertcat(elements{:});
data = [Node, solution];
patch('Faces', elements,'Vertices', data,'FaceColor', 'interp',... 
      'CData', solution / max(abs(solution)));
axis('square')
xlim([min(Node(:,1)) - 0.1, max(Node(:,1)) + 0.1])
ylim([min(Node(:,2)) - 0.1, max(Node(:,2)) + 0.1])
zlim([min(solution) - 0.1, max(solution) + 0.1])
view(60,30)
end