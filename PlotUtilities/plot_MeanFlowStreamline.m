function [x_seed, y_seed] = plot_MeanFlowStreamline(ax, grid_vis, U_vis, V_vis)

X_vis = grid_vis.X;
Y_vis = grid_vis.Y;

hold on

% x_seed = linspace(0.1*min(X_vis(:)), 0.75*max(X_vis(:)), 4);
% y_seed = linspace(min(Y_vis(:)), max(Y_vis(:)), 4);

x_seed = linspace(0.1*min(X_vis(:)), 0.9*max(X_vis(:)), 6);
y_seed = linspace(min(Y_vis(:)), max(Y_vis(:)), 6);


[x_seedmesh, y_seedmesh] = meshgrid(x_seed, y_seed);

grey_vector = [0.15,0.15,0.15];
white_vector = [0.85,0.85,0.85];
cyan_vector = [0,1,1];
h = streamline(ax, X_vis',Y_vis',U_vis',V_vis',x_seedmesh(:),y_seedmesh(:));
set(h,'Color',cyan_vector, 'LineStyle','--', 'LineWidth',2)
h = streamline(ax, X_vis',Y_vis',-U_vis',-V_vis',x_seedmesh(:),y_seedmesh(:));
set(h,'Color',cyan_vector, 'LineStyle','--', 'LineWidth',2)
end