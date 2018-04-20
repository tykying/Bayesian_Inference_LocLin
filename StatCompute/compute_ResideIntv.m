function ResideIntv_Stat = compute_ResideIntv(TransitStat_Cell, Mesh_Struct)

NCells = length(Mesh_Struct.Mesh);
ResideIntv_Stat = zeros(NCells, 4);
for cell_k = 1:NCells
    Transit_Cell = TransitStat_Cell(cell_k);
    
    ResideIntv_Data = [Transit_Cell.ResideIntv(:); Transit_Cell.ResideIntv_I(:); Transit_Cell.ResideIntv_F(:)];
    
    % Columns: Max, Mean, Median, Mode (in alphabetic order)
    ResideIntv_Stat(cell_k, 1) = max(ResideIntv_Data);
    ResideIntv_Stat(cell_k, 2) = mean(ResideIntv_Data);
    ResideIntv_Stat(cell_k, 3) = median(ResideIntv_Data);
    ResideIntv_Stat(cell_k, 4) = mode(ResideIntv_Data);
end

end