function grid_vis = setup_grid_vis(Mesh, RectMesh_Param)

cellij_to_k = RectMesh_Param.cellij_to_k;
Nx_cell = RectMesh_Param.Nx_cell;
Ny_cell = RectMesh_Param.Ny_cell;

X = zeros(Nx_cell, Ny_cell);
Y = zeros(Nx_cell, Ny_cell);
cell_k = zeros(Nx_cell, Ny_cell);
cell_centre = zeros(Nx_cell, Ny_cell, 2);


for jy=1:Ny_cell
    for ix=1:Nx_cell
        cell_k_ij = cellij_to_k(ix, jy);
        cell_centre_ij = Mesh(cell_k_ij).cell_centre;
        
        X(ix, jy) = cell_centre_ij(1);
        Y(ix, jy) = cell_centre_ij(2);
        cell_k(ix, jy) = cell_k_ij;
        cell_centre(ix, jy, :) = cell_centre_ij;
    end
end

grid_vis = struct('X',X, 'Y',Y, 'cell_centre', cell_centre, 'cell_k', cell_k);

end
