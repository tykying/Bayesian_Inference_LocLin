function [TrajPosition_c, Mesh_Struct_c] = label_CoarseGrid(TrajData, RectMesh_Param, Coarsen_factor)
    SpDis_ranges_min = RectMesh_Param.SpDis_ranges_min;
    SpDis_ranges_max = RectMesh_Param.SpDis_ranges_max;

    % Coarsening
    Nx_cell_c = RectMesh_Param.Nx_cell/Coarsen_factor;
    Ny_cell_c = RectMesh_Param.Ny_cell/Coarsen_factor;
    
    assert(mod(Nx_cell_c, 1) == 0);
    assert(mod(Ny_cell_c, 1) == 0);

    
    [Mesh_c, RectMesh_Param_c] = setup_RectMeshStruct(Nx_cell_c, Ny_cell_c, SpDis_ranges_min, SpDis_ranges_max);
    Mesh_Struct_c = struct('Mesh', Mesh_c, 'RectMesh_Param', RectMesh_Param_c);


    [GridInd_c, alphax_c, alphay_c] = positioning_RectMesh(TrajData.x, TrajData.y, RectMesh_Param_c);

    TrajPosition_c = struct('GridInd', GridInd_c, ...
                       'alphax', alphax_c, 'alphay', alphay_c );
end