%% Setup the function space to sample from
% Specifications:

% % Input:
% Nx_cell, Ny_cell: Number of cells in each direction
% SpDis_ranges_min, SpDis_ranges_max: Range over which 

% % Output:  
% (Variables)
% Mesh: a struct that contains info of all the cells, already in cell(k) indexing
% grid_centre_x, grid_centre_y: struct array that contains all individual jumps
% Mesh: a struct that contains info of all the cells, already in cell(k) indexing

% (Function handle)
% cellij_to_k: Map from cell(i,j) to cell(k)
% cellk_to_ij: Inverse Map from cell(i,j) to cell(k)

% Reminder:
% Two cell indexing systems: cell(i,j) and cell(k) => I chose the latter one for generality

%% Computation
function [Mesh, RectMesh_Param] = setup_RectMeshStruct(Nx_cell, Ny_cell, SpDis_ranges_min, SpDis_ranges_max)

Nx_grid = Nx_cell+1;  % number of grid points in x direction
Ny_grid = Ny_cell+1;  % number of grid points in y direction

Ncell_SpDis = Nx_cell*Ny_cell;  % Ncell_SpDis: total number of cells

grid_x = linspace(SpDis_ranges_min(1), SpDis_ranges_max(1), Nx_grid);  % grid points (row vector)
grid_y = linspace(SpDis_ranges_min(2), SpDis_ranges_max(2), Ny_grid);  % grid points

dx = diff(grid_x);
dy = diff(grid_y);

grid_centre_x = 0.5*(grid_x(1:end-1) + grid_x(2:end));
grid_centre_y = 0.5*(grid_y(1:end-1) + grid_y(2:end));

% Quite tricky: Want the reminder to range from 1 to Nx_cell
% Explanation:
%    k = (j-1)*Nx + i <=> k-1 = (j-1)*Nx + (i-1)
% => (k-1) = (i-1) mod Nx, where i-1 ranges from 0 to Nx-1
% => i ranges from 1 to Nx


% Set-up Mesh struct
Mesh = struct('Ind', {}, 'cell_centre', {}, 'cell_dx', {});
Mesh(Ncell_SpDis).cell_centre = zeros(2, 1);
Mesh(Ncell_SpDis).cell_dx = zeros(2, 1);
Mesh(Ncell_SpDis).Ind = zeros(1, 1);

% Indexing: Keep const j, loop over i first.
cell = 0;
for j_cell = 1:Ny_cell
    for i_cell = 1:Nx_cell
        cell = cell + 1;

        Mesh(cell).Ind = cell;
        Mesh(cell).cell_centre = [grid_centre_x(i_cell); grid_centre_y(j_cell)];
        Mesh(cell).cell_dx = [dx(i_cell); dy(j_cell)];
        % Mesh(cell).connectivity   %% Do not need it for the moment; Perhaps needed for adaptivity
    end
end


% Unused functions: For future reference
% Two cell indexing systems: cell(i,j) and cell(k)
% Set-up two function that converts (i,j) into k
cellij_to_k = @(i, j) (j-1)*Nx_cell + i;
cellk_to_ij = @(k) [mod((k-1), Nx_cell)+1, floor((k-1)/Nx_cell)+1];  

bin_sizes = (SpDis_ranges_max - SpDis_ranges_min)./([Nx_cell, Ny_cell]);

RectMesh_Param = struct( 'Nx_cell', Nx_cell, 'Ny_cell', Ny_cell, ...
                   'SpDis_ranges_min', SpDis_ranges_min, ...
                   'SpDis_ranges_max', SpDis_ranges_max, ...
                   'grid_x', grid_x, 'grid_y', grid_y, ...
                   'bin_sizes', bin_sizes, ...
                   'cellij_to_k', cellij_to_k,  'cellk_to_ij', cellk_to_ij);

end