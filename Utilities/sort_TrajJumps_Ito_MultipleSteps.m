% Associate the jumps to the starting points
% Used for exact likelihood from EA
function [TrajJumps, NJumps_cell] = sort_TrajJumps_Ito_MultipleSteps(TrajData, TrajPosition, Mesh, Nlevel)
    NCells = length(Mesh);
    
    % A (Determine Blackbox) scheme to take the multiple steps
    % Do not want to double use the data
    % Do not want to choose correlated data
    % Want evening sampling in the entire time-interval
    
    % Nlevel: should not exceeds 5 for efficent use of observation data 
    % Nlevel: 4 or 5 
    MaxPow = (Nlevel-1);
    MaxSep = 2^MaxPow;
    
    % A (Determine Blackbox) scheme to take the multiple steps
    % Do not want to double use the data
    % Do not want to choose correlated data
    % Want evening sampling in the entire time-interval
    level_endpts = obtain_level_endpts(Nlevel);
    
    %assert(length(unique(level_endpts(:))) == size(level_endpts, 1) * size(level_endpts, 2));
    
    % Verify the raw sampling interval is uniform
    h_list = diff(TrajData.ts_list);
    assert(range(h_list) < 0.1);
    
    [Nts, nparticles] = size(TrajData.x);
    
    tIndref_filter = 1:MaxSep:(Nts-MaxSep);
    for level = 1:Nlevel
        tInd0_filter = tIndref_filter + level_endpts(level, 1);
        tInd1_filter = tIndref_filter + level_endpts(level, 2);
        
        h_filter = TrajData.ts_list(tInd1_filter)-TrajData.ts_list(tInd0_filter);
        [X, h] = meshgrid(ones(1, nparticles), h_filter);  % No longer matter whether h_filter is a row or column matrix       
        
        diffx = TrajData.x(tInd1_filter, :) - TrajData.x(tInd0_filter, :);
        diffy = TrajData.y(tInd1_filter, :) - TrajData.y(tInd0_filter, :);
   
        x0 = TrajData.x(tInd0_filter, :);
        y0 = TrajData.y(tInd0_filter, :);
        
        alphax = TrajPosition.alphax(tInd0_filter, :);
        alphay = TrajPosition.alphay(tInd0_filter, :);
        
        GridInd = TrajPosition.GridInd(tInd0_filter, :);
        
        for cell_k = 1:NCells
            if (mod(cell_k, 8*4) == 0)
                disp(['Level ', num2str(level), ': Sorting Jumps in cell', num2str(cell_k)]);
            end
            
            CellCentre = Mesh(cell_k).cell_centre;
            GridInd_CellData = (GridInd == cell_k);
            
            % Shift X_0 relative to the cell_centre
            x0_CellData = x0(GridInd_CellData) - CellCentre(1);
            y0_CellData = y0(GridInd_CellData) - CellCentre(2);
            
            alphax_CellData = alphax(GridInd_CellData);
            alphay_CellData = alphay(GridInd_CellData);
            
            diffx_CellData = diffx(GridInd_CellData);
            diffy_CellData = diffy(GridInd_CellData);
            
            h_CellData = h(GridInd_CellData);

            % Make sure a column vector is output
            if size(x0_CellData, 2) ~= 1
                x0_CellData = x0_CellData';
                y0_CellData = y0_CellData';
                alphax_CellData = alphax_CellData';
                alphay_CellData = alphay_CellData';
                diffx_CellData = diffx_CellData';
                diffy_CellData = diffy_CellData';
                h_CellData = h_CellData';
            end
            
            TrajJumps(cell_k, level) = struct('x0', x0_CellData, ...
                'y0', y0_CellData, ...
                'alphax', alphax_CellData, ...
                'alphay', alphay_CellData, ...
                'diffx', diffx_CellData, ...
                'diffy', diffy_CellData, ...
                'h', h_CellData);
        end
    end
    
    NJumps_cell = zeros(NCells, Nlevel);
    for level = 1:Nlevel
        for cell_k = 1:NCells
            NJumps_cell(cell_k, level) = NJumps_cell(cell_k, level) + length(TrajJumps(cell_k, level).x0);
        end
    end
    
    
    % Deterministic Blackbox to take the multiple steps
    function level_endpts = obtain_level_endpts(Nlevel)
        MaxPow = (Nlevel-1);
        MaxSep = 2^MaxPow;
        
        level_endpts = zeros(Nlevel, 2);
        level_refpts = zeros(Nlevel, 1);
        
        level_refpts(1) = MaxSep;
        level_refpts_prev = level_refpts(1);
        diffFac = 2^(MaxPow-2)*2;
        for level = 2:(round(MaxPow/2)+1)
            diffFac = diffFac/2;
            level_refpts(level) = level_refpts_prev - diffFac;
            level_refpts_prev = level_refpts(level);
        end
        if mod(MaxPow, 2) == 0
            diffFac = diffFac;
        else
            diffFac = diffFac*2;
        end
        for level = (round(MaxPow/2)+2):Nlevel
            level_refpts(level) = level_refpts_prev - diffFac-1;
            level_refpts_prev = level_refpts(level);
            diffFac = diffFac*2;
        end
        
        level_width = MaxSep;
        level_endpts(1,:) = [0, MaxSep];
        for level = 2:(round(MaxPow/2)+1)
            if mod(level, 2) == 0
                % level_refpts as ending point:
                level_width = level_width/2;
                Idx_f = level_refpts(level);
                level_endpts(level, :) = [Idx_f-level_width, Idx_f];
            elseif mod(level, 2) == 1
                % When choosing starting point:
                % Middle of last starting point (2 levels before) and mid-point of the domain
                level_width = level_width/2;
                Idx_i = level_refpts(level);
                level_endpts(level, :) = [Idx_i, Idx_i+level_width];
            end
        end
        for level = (round(MaxPow/2)+2):Nlevel
            if mod(level, 2) == 0
                % level_refpts as ending point:
                level_width = level_width/2;
                Idx_f = level_refpts(level)+2;
                level_endpts(level, :) = [Idx_f-level_width, Idx_f];
            elseif mod(level, 2) == 1
                % When choosing starting point:
                % Middle of last starting point (2 levels before) and mid-point of the domain
                level_width = level_width/2;
                Idx_i = level_refpts(level);
                level_endpts(level, :) = [Idx_i, Idx_i+level_width];
            end
        end
        
        if Nlevel == 2
            disp('Special case for Nlevel = 2')
            level_endpts = [0 2; 1, 2];
        end
    end
end