function [GridInd, alphax, alphay] = positioning_RectMesh(x,y, RectMesh_Param)
    cellij_to_k = RectMesh_Param.cellij_to_k;
    bin_sizes = RectMesh_Param.bin_sizes;
    SpDis_ranges_min = RectMesh_Param.SpDis_ranges_min;
    Nx_cell = RectMesh_Param.Nx_cell;
    Ny_cell = RectMesh_Param.Ny_cell;
    
    [Nts, nparticles] = size(x);
    
    GridInd = zeros(Nts, nparticles);
    
    alphax = zeros(Nts, nparticles);
    alphay = zeros(Nts, nparticles);

    % Locate the particle in each bin
    for i_t = 1: Nts
        [x_CellInd, x_Cellalpha] = search_sorted_uniform(x(i_t,:), bin_sizes(1), SpDis_ranges_min(1));
        [y_CellInd, y_Cellalpha] = search_sorted_uniform(y(i_t,:), bin_sizes(2), SpDis_ranges_min(2));
        
        % Filter out x_CellInd: all out-of-bound are labelled as -1
        x_CellInd(x_CellInd > Nx_cell) = NaN; x_CellInd(x_CellInd < 1) = NaN;
        y_CellInd(y_CellInd > Ny_cell) = NaN; y_CellInd(y_CellInd < 1) = NaN;
        
        CellIndk_part = cellij_to_k(x_CellInd, y_CellInd);
        
        GridInd(i_t, :) = CellIndk_part;
        
        alphax(i_t, :) = x_Cellalpha;
        alphay(i_t, :) = y_Cellalpha;
    end
end
