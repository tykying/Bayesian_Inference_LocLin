function [Field_f_Itp, X_f, Y_f] = griddata_interpolation_pwl(Field_c_Mean, Field_c_Gradx, Field_c_Grady, RectMesh_Param, FineGrid_Resolution)
% Ensure right dimension of data and RectMesh_Param
assert(RectMesh_Param.Nx_cell == size(Field_c_Mean, 1));
assert(RectMesh_Param.Ny_cell == size(Field_c_Mean, 2));
assert(all(size(Field_c_Mean) == size(Field_c_Gradx)));
assert(all(size(Field_c_Mean) == size(Field_c_Grady)));

% Assert FineGrid_Resolution-1 to be power of 2; 
% In fact NOT necessary
%assert(mod(log2(FineGrid_Resolution-1), 1) == 0);

SpDis_ranges_min = RectMesh_Param.SpDis_ranges_min;
SpDis_ranges_max = RectMesh_Param.SpDis_ranges_max;

% Set up coarse-grid
grid_x_c = RectMesh_Param.grid_x;
grid_y_c = RectMesh_Param.grid_y;

x_c = 0.5*(grid_x_c(1:end-1) + grid_x_c(2:end));
y_c = 0.5*(grid_y_c(1:end-1) + grid_y_c(2:end));

[X_c, Y_c] = meshgrid(x_c, y_c);
X_c = X_c'; Y_c =Y_c';  % So that same orientiation as grid_vis.X;

% Set up fine-grid
x_f = linspace(SpDis_ranges_min(1), SpDis_ranges_max(1), FineGrid_Resolution);
y_f = linspace(SpDis_ranges_min(2), SpDis_ranges_max(2), FineGrid_Resolution);

[X_f, Y_f] = meshgrid(x_f, y_f);
X_f = X_f'; Y_f =Y_f';  % So that same orientiation as grid_vis.X;

% Assoicate a list of GridInd to each point on 1D line
% Enable the point to lie exactly on the Grid
x_f_AssoGridInd = assoicate_points_with_GridInd_uniform(x_f, RectMesh_Param.bin_sizes(1), RectMesh_Param.grid_x(1), RectMesh_Param.Nx_cell);
y_f_AssoGridInd = assoicate_points_with_GridInd_uniform(y_f, RectMesh_Param.bin_sizes(2), RectMesh_Param.grid_y(1), RectMesh_Param.Ny_cell);

% Set-up a 2D Array of List; Each list contains the assoicated GridInd 
AssoGridInd = zeros(length(x_f), length(y_f), size(x_f_AssoGridInd, 1)*size(y_f_AssoGridInd, 1));
for j_f = 1:length(y_f)
    for i_f = 1:length(x_f)
        % At a given point, obtain all the assoicated grids
        [x_f_AGI_ij, y_f_AGI_ij] = meshgrid(x_f_AssoGridInd(:, i_f), y_f_AssoGridInd(:, j_f));
        xy_f_AGI_k = RectMesh_Param.cellij_to_k(x_f_AGI_ij, y_f_AGI_ij);
        AssoGridInd(i_f, j_f, :) = xy_f_AGI_k(:);
    end
end

% For easier interpretation: y = A(x-x_bar) + b
A_x = Field_c_Gradx;
A_y = Field_c_Grady;
b = Field_c_Mean;

Field_f_Itp = zeros(size(AssoGridInd));
% Assume cell centered variable
% y = A(x-x_bar) + b; x_bar: cell centre of the mesh
for j_f = 1:size(Field_f_Itp, 2)
    for i_f = 1:size(Field_f_Itp, 1)
        for GridInd_iterator = 1:size(Field_f_Itp, 3)  % loop over all 4 values
            AssoGridInd_k = AssoGridInd(i_f, j_f, GridInd_iterator);
            ij_c = RectMesh_Param.cellk_to_ij(AssoGridInd_k);
            i_c = ij_c(1); j_c = ij_c(2);
            
            A_ij = [A_x(i_c, j_c), A_y(i_c, j_c)];
            b_ij = b(i_c,j_c);
            
            % x_bar = X_c_ij = Cell centre of coarse grid
            X_c_ij = [X_c(i_c, j_c); Y_c(i_c, j_c)];  % Cell centre of coarse grid
            X_f_ij = [X_f(i_f, j_f); Y_f(i_f, j_f)];  % Coordinates of fine grid point
            % A(x-x_bar)+b
            Y_ij = A_ij*(X_f_ij-X_c_ij) + b_ij;
            
            Field_f_Itp(i_f,j_f, GridInd_iterator) = Y_ij;
        end
    end
end


% Post-process:
% if a fine grid point is only associated to one coarse grid point,
% no need to output all four values
Overlapped_GridPoint = 0;
for j_f = 1:length(y_f)
    for i_f = 1:length(x_f)
        if (length(unique(AssoGridInd(i_f, j_f, :))) ~= 1)
            Overlapped_GridPoint = 1;
        end
    end
end

if Overlapped_GridPoint == 0
    Field_f_Itp = squeeze(Field_f_Itp(:, :, 1));
end

end
%Range_Discont_REL = range(Field_f_Itp, 3)./mean(Field_f_Itp, 3);