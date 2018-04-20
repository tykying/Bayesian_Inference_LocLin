%% Setup Rectangular Structured Grid/Mesh
if (contains(veloc_Profile, 'QGM2')==1)  % =0: Matched
    SpDis_ranges_min = [0, 0]; SpDis_ranges_max = [L, L];
    Nx_cell = Nx_cell_ARG; Ny_cell = Nx_cell_ARG;
else
    [range_min, range_max] = LocateBoundaries(x, y);
    
    if (strcmp(veloc_Profile, 'linear_shear'))  % Ad-hoc
        disp('AD-HOC: set SpDis_L=2 in Script_SetupRectMesh.m');
        SpDis_L = 2;
    elseif (strcmp(veloc_Profile, 'taylor_green_noisy'))
        disp('AD-HOC: set SpDis_L=4 in Script_SetupRectMesh.m');
        SpDis_L = 4;
    elseif (strcmp(veloc_Profile, 'childress_soward'))
        disp('AD-HOC: set SpDis_L=4 in Script_SetupRectMesh.m');
        SpDis_L = 4;
    elseif (strcmp(veloc_Profile, 'Mixed'))
        disp('AD-HOC: set SpDis_L=4 in Script_SetupRectMesh.m');
        SpDis_L = 4;    
    else
        % Consider only square domain
        SpDis_L = max([abs(range_min), abs(range_max)]);
    end
    SpDis_ranges_min = [-SpDis_L, -SpDis_L];
    SpDis_ranges_max = [SpDis_L, SpDis_L];
    Nx_cell = Nx_cell_ARG; Ny_cell = Nx_cell_ARG;
    
    if (contains(veloc_Profile, 'cosine'))
        disp('AD-HOC: set domain in Script_SetupRectMesh.m');
        SpDis_L = round(max([abs(range_min), abs(range_max)]));
        
        SpDis_ranges_min = [-2^round(log2(SpDis_L)+0.5) -1];
        SpDis_ranges_max = [2^round(log2(SpDis_L)+0.5), 1];
        
        Nx_cell = 1;
        Ny_cell = Nx_cell_ARG;
    end
end

[Mesh, RectMesh_Param] = setup_RectMeshStruct(Nx_cell, Ny_cell, SpDis_ranges_min, SpDis_ranges_max);
Mesh_Struct = struct('Mesh', Mesh, 'RectMesh_Param', RectMesh_Param);

bin_sizes = RectMesh_Param.bin_sizes;