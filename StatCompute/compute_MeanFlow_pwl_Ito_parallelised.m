function [MeanFlow, EffectiveNJumps, EddiesData, MeanFlowData] = compute_MeanFlow_pwl_Ito_parallelised(TrajData, TrajPosition, Mesh_Struct)
type = 'pwl';

Mesh = Mesh_Struct.Mesh;
NCells = length(Mesh);

x = TrajData.x;
y = TrajData.y;
ts_list = TrajData.ts_list;

GridInd = TrajPosition.GridInd;
GridInd(end, :) = [];  % Ito: shift from t_n to t_n; Remove the last row for displacement

diffx = diff(x, 1);
diffy = diff(y, 1);

difft = zeros(size(diffx));
for part = 1:size(diffx, 2)
    difft(:, part) = diff(ts_list, 1);
end

% Partition u, v into 16 equal pieces
NNodes = 4;
Node_PartList = round(linspace(1, size(diffx, 2), NNodes+1));
for node = 1:NNodes
    PartInd_i = Node_PartList(node);
    PartInd_f = Node_PartList(node+1);
    
    x_node = x(:, PartInd_i:PartInd_f);
    y_node = y(:, PartInd_i:PartInd_f);
    
    difft_node = difft(:, PartInd_i:PartInd_f);
    
    diffx_node = diffx(:, PartInd_i:PartInd_f);
    diffy_node = diffy(:, PartInd_i:PartInd_f);
    GridInd_node = GridInd(:, PartInd_i:PartInd_f);
    
    Node_TrajJumps(node) = struct('x', x_node, 'y', y_node, ...
        'diffx', diffx_node, 'diffy', diffy_node, ...
        'difft', difft_node, 'GridInd', GridInd_node);
end
clear var x y diffx diffy GridInd 
clear var x_node y_node u_node v_node difft_node GridInd_node

% Parallelise: Sorting into cells
Node_TrajJumps_inCell(NCells) = struct('x', [], 'y', [], ...
        'diffx', [], 'diffy', [], 'difft', [], 'GridInd', []);
    
% Nested Structure
for node = 1:NNodes
    NodeCell_TrajJumps{node} = Node_TrajJumps_inCell;
end

parfor node = 1:NNodes
    TrajJumps = Node_TrajJumps(node);
    x = TrajJumps.x;
    y = TrajJumps.y;
    diffx = TrajJumps.diffx;
    diffy = TrajJumps.diffy;
    difft = TrajJumps.difft;
    GridInd = TrajJumps.GridInd;
    
    Node_TrajJumps_inCell(NCells) = struct('x', [], 'y', [], ...
        'diffx', [], 'diffy', [], 'difft', [], 'GridInd', []);
    
    for cell_k = 1:NCells  
        if (mod(log2(cell_k), 1) == 0)
            fprintf('Parallel sorting in cell %d \n', cell_k);
        end
        
        Node_TrajJumps_inCell(cell_k) = sort_JumpsPosition_inCell(x, y, diffx, diffy, GridInd, difft, cell_k);
    end
    
    NodeCell_TrajJumps{node} = Node_TrajJumps_inCell;
end
clear var Node_TrajJumps;
disp('Finished sorting');

TrajJumps(NCells) = struct('x', [], 'y', [], ...
        'diffx', [], 'diffy', [], 'difft', [], 'GridInd', []);
% Sequentially Merge all nodes into cells
for cell_k = 1:NCells
    if (mod(log2(cell_k), 1) == 0)
        fprintf('Sequential Reording in cell %d \n', cell_k);
    end
    
    diffx = [];
    diffy = [];
    x = [];
    y = [];
    GridInd = [];
    difft = [];
    
    for node = 1:NNodes
        diffx = [diffx; NodeCell_TrajJumps(node, cell_k).diffx];
        diffy = [diffy; NodeCell_TrajJumps(node, cell_k).diffy];
        difft = [difft; NodeCell_TrajJumps(node, cell_k).difft];

        x = [x; NodeCell_TrajJumps(node, cell_k).x];
        y = [y; NodeCell_TrajJumps(node, cell_k).y];
        GridInd = [GridInd; NodeCell_TrajJumps(node, cell_k).GridInd];
    end
    
    TrajJumps(cell_k).x = x;
    TrajJumps(cell_k).y = y;
    TrajJumps(cell_k).u = diffx;
    TrajJumps(cell_k).v = diffy;
    TrajJumps(cell_k).GridInd = GridInd;
end
clear var u v x y GridInd

% % Parallisation by partition of data
% for node = 1:NNodes
%     % Ito: Place the 'velocities' at the initial positions
%     % Estimated velocity fields
%     U_m = zeros(NCells, 2);
%     GradU_m = zeros(NCells, 2, 2);  % First 2: 2 components; Second 2: u and v
%     
%     % Store the decomposed velocity data
%     u_mean = zeros(size(u));
%     v_mean = zeros(size(v));
%     u_eddy = zeros(size(u));
%     v_eddy = zeros(size(v));
%     
%     EffectiveNJumps = zeros(NCells, 1);
%     
%     for cell_k = 1:NCells
%         if mod(log2(cell_k), 2) == 0
%             fprintf('Fitting in cell %i. \n', cell_k);
%         end
%         
%         GridInd_inCell = (GridInd == cell_k);
%         
%         CellCen = Mesh(cell_k).cell_centre';
%         
%         %   %%% Only when spatial structure is involved
%         X_inCell = [x(GridInd_inCell), y(GridInd_inCell)];
%         
%         u_inCell = u(GridInd_inCell);
%         v_inCell = v(GridInd_inCell);
%         
%         % Decomposing instantaneous velocity by linear fitting
%         [u_sf, u_goodness] = fit(X_inCell, u_inCell, 'poly11');
%         [v_sf, v_goodness] = fit(X_inCell, v_inCell, 'poly11');
%         
%         u_mean(GridInd_inCell) = feval(u_sf, X_inCell);
%         v_mean(GridInd_inCell) = feval(v_sf, X_inCell);
%         
%         u_eddy(GridInd_inCell) = u_inCell - u_mean(GridInd_inCell);
%         v_eddy(GridInd_inCell) = v_inCell - v_mean(GridInd_inCell);
%         
%         % Obtain cell-averaged values
%         U_m(cell_k, :) = [feval(u_sf, CellCen), feval(v_sf, CellCen)];
%         
%         GradU_m(cell_k, :, 1) = [u_sf.p10, u_sf.p01];
%         GradU_m(cell_k, :, 2) = [v_sf.p10, v_sf.p01];
%         
%         EffectiveNJumps(cell_k) = length(u_inCell);
%     end
% end
% 
% % Idea for future remark:
% % Parallel and orthogonal to mean flow direction?
% 
% % For Output
% % MeanFlowData = struct('u', u_mean, 'v', v_mean, 'ds', ds_mean);
% % EddiesData = struct('u', u_eddy, 'v', v_eddy, 'ds', ds_eddy);
% MeanFlowData = struct('u', u_mean, 'v', v_mean);
% EddiesData = struct('u', u_eddy, 'v', v_eddy);
% 
% MeanFlow = struct('U_m', U_m, 'GradU_m', GradU_m, 'type', type);
end
