function [ Data_Cell_k ] = convert_cellij_to_cellk( Data_Cell_ij, RectMesh_Param )

Nx_cell = RectMesh_Param.Nx_cell;
Ny_cell = RectMesh_Param.Ny_cell;

[Nx_cell, Ny_cell, NDataCate] = size(Data_Cell_ij);


NCells = Nx_cell*Ny_cell;

Data_Cell_k = zeros(NCells, NDataCate);

for DataCate = 1 : NDataCate
    for cell_i = 1 : Nx_cell
        for cell_j = 1 : Ny_cell
        cell_k = RectMesh_Param.cellij_to_k(cell_i, cell_j);
        
        Data_Cell_k(cell_k, DataCate) = Data_Cell_ij(cell_i, cell_j, DataCate);
        end
    end
end

% Why did I ever do this?? Keeping the same structure should be better?
% Data_Cell_ij = squeeze(Data_Cell_ij);

end