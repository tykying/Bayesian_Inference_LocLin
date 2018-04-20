function [x_seed, y_seed] = plot_VelFieldStreamline(veloc_fldStruct, grid_vis)
X_vis = grid_vis.X;
Y_vis = grid_vis.Y;

U_vis = zeros(size(X_vis));
V_vis = zeros(size(X_vis));

for i = 1:size(X_vis, 1)
    for j = 1:size(X_vis, 2)
        u = veloc_fldStruct.u(X_vis(i,j), Y_vis(i,j));
        U_vis(i,j) = u(1);
        V_vis(i,j) = u(2);
    end
end

hold on

% x_seed = linspace(0.1*min(X_vis(:)), 0.75*max(X_vis(:)), 4);
% y_seed = linspace(min(Y_vis(:)), max(Y_vis(:)), 4);

x_seed = linspace(0.1*min(X_vis(:)), 0.75*max(X_vis(:)), 5);
y_seed = linspace(min(Y_vis(:)), max(Y_vis(:)), 5);


x_seed = linspace(min(X_vis(:)), max(X_vis(:)), 10);
y_seed = linspace(min(Y_vis(:)), max(Y_vis(:)), 10);

[x_seedmesh, y_seedmesh] = meshgrid(x_seed, y_seed);

grey_vector = [0.2,0.2,0.2];
white_vector = [0.85,0.85,0.85];
h = streamline(X_vis',Y_vis',U_vis',V_vis',x_seedmesh(:),y_seedmesh(:));
set(h,'Color',grey_vector, 'LineStyle','--')
%h = streamline(X_vis',Y_vis',-U_vis',-V_vis',x_seedmesh(:),y_seedmesh(:));
%set(h,'Color',grey_vector, 'LineStyle','--')

end