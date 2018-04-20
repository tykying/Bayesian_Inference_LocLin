function theta2D = thetaStruct_to_Rect2DAarray(theta, RectMesh_Param)

% Should check theta Struct is of the same form
%theta_sum = struct('U',{},'V',{},'K',{},'GradU',{},'GradV',{},'GradK',{});
Nx_cell = RectMesh_Param.Nx_cell;
Ny_cell = RectMesh_Param.Ny_cell;
cellij_to_k = RectMesh_Param.cellij_to_k;
cellk_to_ij = RectMesh_Param.cellk_to_ij;

thetaFields = fieldnames(theta(1));
N_Fields = numel(thetaFields);

% Create another object theta with same fields
theta2D = theta;

for j = 1: N_Fields
    FldName = thetaFields{j};
    
    [RawDataNCell, RawDataDim] = size(theta(1).(FldName));  % First dimension = number of cells
    theta2D.(FldName) = zeros(Nx_cell, Ny_cell, RawDataDim);
end

for j = 1: N_Fields
    FldName = thetaFields{j};
    for cell_j = 1 : Ny_cell
        for cell_i = 1 : Nx_cell
            cell_k = cellij_to_k(cell_i, cell_j);
        
            theta2D.(FldName)(cell_i, cell_j, :) = theta.(FldName)(cell_k, :);
        end
    end
end


end