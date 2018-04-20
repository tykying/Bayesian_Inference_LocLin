% Reshape DATA of size=(cell_k, NDataCate) into (cell_i, cell_j, NDataCate)
function [ Data_Cell_ij ] = convert_cellk_to_cellij( Data_Cell_k, RectMesh_Param )

Nx_cell = RectMesh_Param.Nx_cell;
Ny_cell = RectMesh_Param.Ny_cell;

[NCells, NDataCate] = size(Data_Cell_k);


assert(NCells == Nx_cell*Ny_cell);

Data_Cell_ij = zeros(Nx_cell, Ny_cell, NDataCate);

for DataCate = 1 : NDataCate
    for cell_k = 1 : NCells
        cell_ij = RectMesh_Param.cellk_to_ij(cell_k);
        cell_i = cell_ij(1);
        cell_j = cell_ij(2);
        
        Data_Cell_ij(cell_i, cell_j, DataCate) = Data_Cell_k(cell_k, DataCate);
    end
end

% Why did I ever do this?? Keeping the same structure should be better?
% Data_Cell_ij = squeeze(Data_Cell_ij);

end