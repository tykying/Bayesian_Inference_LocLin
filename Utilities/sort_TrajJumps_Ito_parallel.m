% Associate the jumps to the starting points
% Used for exact likelihood from EA


%% PARALLEL VERSION DOES NOT REALLY WORK
% Fighting for shared data

function [TrajJumps] = sort_TrajJumps_Ito_parallel(TrajData, TrajPosition, Mesh)
    NCells = length(Mesh);
    
    % Assumed uniform sampling interval  
    diffx_Global = diff(TrajData.x, 1);
    diffy_Global = diff(TrajData.y, 1);

    h_list = diff(TrajData.ts_list);
    if (size(h_list, 2) ~= 1)  % QGM2 and Matlab: I store them in different order
        h_list = h_list';
    end
    assert(size(h_list, 2) == 1);
    h_Global = (h_list)*ones(1, size(TrajData.x, 2));  % Outer product: generate matrix h
    
    % Ito velocity: Filter out the last time-step
    x0_Global = TrajData.x(1:end-1, :);
    y0_Global = TrajData.y(1:end-1, :);
    
    GridInd_Global = TrajPosition.GridInd(1:end-1, :);
    alphax_Global = TrajPosition.alphax(1:end-1, :);
    alphay_Global = TrajPosition.alphay(1:end-1, :);
    
    assert(all(size(GridInd_Global) == size(diffx_Global)));
    
    %TrajJumps_sorted = {};
    %for cell_k = 1:NCells
    parfor cell_k = 1:NCells
        % Obtain local copies [Appears to be extremely slow]
        diffx = diffx_Global;
        diffy = diffy_Global;
        
        x0 = x0_Global;
        y0 = y0_Global;
        
        h = h_Global;
        
        GridInd = GridInd_Global;
        alphax = alphax_Global;
        alphay = alphay_Global;
        
        if (mod(cell_k, 8*4) == 0)
            disp(['Sorting Jumps in cell', num2str(cell_k)]);
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
        end
        if size(h_CellData, 2) ~= 1
            h_CellData = h_CellData';
        end
        
        TrajJumps_Local = struct('x0', x0_CellData, ...
            'y0', y0_CellData, ...
            'alphax', alphax_CellData, ...
            'alphay', alphay_CellData, ...
            'diffx', diffx_CellData, ...
            'diffy', diffy_CellData, ...
            'h', h_CellData);
        
%         TrajJumps(cell_k) = struct('x0', x0_CellData, ...
%             'y0', y0_CellData, ...
%             'alphax', alphax_CellData, ...
%             'alphay', alphay_CellData, ...
%             'diffx', diffx_CellData, ...
%             'diffy', diffy_CellData, ...
%             'h', h_CellData);
    end
    
    
end