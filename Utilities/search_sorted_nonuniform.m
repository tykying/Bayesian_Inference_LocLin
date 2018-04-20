function [x_CellInd, x_alpha] = search_sorted_nonuniform(x, grid_boundary)
%% Locate in which grid x lies in, and return the relative position in the grid

% Specifications:

% Input: (all scalars)
% x: one coordinate of a particle
% grid_boundary: boundary coordinates of non-uniform grid, assumed to be sorted in increasing order

% Output:
% x_CellInd: tell in which cell i, x lies in
% x_alpha_cell: coresponding relative position in the cell

%% computation
dx = diff(grid_boundary);

x_CellInd = find(grid_boundary > x, 1) - 1;
x_alpha = (x - grid_boundary(x_CellInd))./dx(x_CellInd);

end
