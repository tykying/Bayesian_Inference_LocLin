function [MeanFlow, EffectiveNJumps, EddiesData, MeanFlowData] = compute_MeanFlow_pwl_Ito(TrajData, TrajPosition, Mesh_Struct)
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
difft = diff(ts_list, 1);

u = diffx./difft;
v = diffy./difft;
% ds_full = sqrt(u.^2 + v.^2);

% Ito: Place the 'velocities' at the initial positions
% Estimated velocity fields
U_m = zeros(NCells, 2);
GradU_m = zeros(NCells, 2, 2);  % First 2: 2 components; Second 2: u and v

% Store the decomposed velocity data
u_mean = zeros(size(u));
v_mean = zeros(size(v));
u_eddy = zeros(size(u));
v_eddy = zeros(size(v));

EffectiveNJumps = zeros(NCells, 1);


% Partition u, v into 16 equal pieces
NNodes = 16;
for node = 1:NNodes
    
    TrajJumps_node(node) = struct('u', u, 'v', v, ...
        'difft', difft, 'GridInd', GridInd);
end
clear var diffx diffy

% %%% Only when spatial structure is involved
% x(end, :) = [];
% y(end, :) = [];
for cell_k = 1:NCells
    if mod(log2(cell_k), 2) == 0
        fprintf('Fitting in cell %i. \n', cell_k);
    end
    
    GridInd_inCell = (GridInd == cell_k);
    
    CellCen = Mesh(cell_k).cell_centre';
    
%   %%% Only when spatial structure is involved
    X_inCell = [x(GridInd_inCell), y(GridInd_inCell)];
    
    u_inCell = u(GridInd_inCell);
    v_inCell = v(GridInd_inCell);
    
    % Decomposing instantaneous velocity by linear fitting
    [u_sf, u_goodness] = fit(X_inCell, u_inCell, 'poly11');
    [v_sf, v_goodness] = fit(X_inCell, v_inCell, 'poly11');
    
    u_mean(GridInd_inCell) = feval(u_sf, X_inCell);
    v_mean(GridInd_inCell) = feval(v_sf, X_inCell);

    u_eddy(GridInd_inCell) = u_inCell - u_mean(GridInd_inCell);
    v_eddy(GridInd_inCell) = v_inCell - v_mean(GridInd_inCell);
    
    % Obtain cell-averaged values
    U_m(cell_k, :) = [feval(u_sf, CellCen), feval(v_sf, CellCen)];

    GradU_m(cell_k, :, 1) = [u_sf.p10, u_sf.p01];
	GradU_m(cell_k, :, 2) = [v_sf.p10, v_sf.p01];
    
    EffectiveNJumps(cell_k) = length(u_inCell);
end

% Idea for future remark:
% Parallel and orthogonal to mean flow direction?

% For Output
% MeanFlowData = struct('u', u_mean, 'v', v_mean, 'ds', ds_mean);
% EddiesData = struct('u', u_eddy, 'v', v_eddy, 'ds', ds_eddy);
MeanFlowData = struct('u', u_mean, 'v', v_mean);
EddiesData = struct('u', u_eddy, 'v', v_eddy);

MeanFlow = struct('U_m', U_m, 'GradU_m', GradU_m, 'type', type);
end
