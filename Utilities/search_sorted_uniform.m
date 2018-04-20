function [x_CellInd, x_alpha] = search_sorted_uniform(x, bin_size, grid_startpt)
%% Locate in which grid x lies in, and return the relative position in the grid

% Specifications:

% Input: (all scalars)
% x: Normalised Position (relative to the domain length)
% grid_boundary: boundary coordinates of uniform grid

% Output:
% x_CellInd: tell in which cell i, x lies in
% x_alpha: coresponding relative position in the cell

%% computation

x = x-grid_startpt;
x(x<0) = nan;  % To prevent bug
x_CellInd = ceil(x/bin_size);  % Equivalent to floor(x/bin_size) + 1;
x_alpha = (x - (x_CellInd-1)*bin_size)/bin_size;

end
