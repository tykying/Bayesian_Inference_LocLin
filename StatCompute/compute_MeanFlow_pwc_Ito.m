function [MeanFlow, EffectiveNJumps, EddiesData, MeanFlowData] = compute_MeanFlow_pwc_Ito(TrajData_Positioned, Mesh_Struct)
type = 'pwc';

NCells = length(Mesh_Struct.Mesh);

x = TrajData_Positioned.x;
y = TrajData_Positioned.y;
ts_list = TrajData_Positioned.ts_list;
GridInd = TrajData_Positioned.GridInd;

diffx = diff(x, 1);
diffy = diff(y, 1);
difft = diff(ts_list, 1);

u = diffx./difft';
v = diffy./difft';
% ds_full = sqrt(u.^2 + v.^2);

% Ito: Place the 'velocities' at the initial positions
u_m = zeros(NCells, 1);
v_m = zeros(NCells, 1);
U_m = zeros(NCells, 2);
GradU_m = zeros(NCells, 2, 2);  % First 2: 2 components; Second 2: u and v

% Store the decomposed velocity data
u_mean = zeros(size(u));
v_mean = zeros(size(v));
u_eddy = zeros(size(u));
v_eddy = zeros(size(v));

EffectiveNJumps = zeros(NCells, 1);

% %%% Only when spatial structure is involved
% x(end, :) = [];
% y(end, :) = [];
for cell_k = 1:NCells
    GridInd_inCell = (GridInd == cell_k);
    GridInd_inCell(end, :) = [];
    
%   %%% Only when spatial structure is involved
%    x_inCell = x(GridInd_inCell);
%    y_inCell = y(GridInd_inCell);
    
    u_inCell = u(GridInd_inCell);
    v_inCell = v(GridInd_inCell);
    
    u_m(cell_k) = mean(u_inCell);
    v_m(cell_k) = mean(v_inCell);
    
    % Decomposing instantaneous velocity
    u_mean(GridInd_inCell) = u_m(cell_k);
    v_mean(GridInd_inCell) = v_m(cell_k);
    
    u_eddy(GridInd_inCell) = u_inCell - u_m(cell_k);
    v_eddy(GridInd_inCell) = v_inCell - v_m(cell_k);
    
    EffectiveNJumps(cell_k) = length(u_inCell);
end

% Compute other eddies velocities components
ds_mean = sqrt(u_mean.^2 + v_mean.^2);
ds_eddy = sqrt(u_eddy.^2 + v_eddy.^2);
% Parallel and orthogonal to mean flow direction?


% U_m: average value in the grid point
U_m = [u_m, v_m];

% For Output
MeanFlowData = struct('u_mean', u_mean, 'v_mean', v_mean, 'ds_mean', ds_mean);
EddiesData = struct('u_eddy', u_eddy, 'v_eddy', v_eddy, 'ds_eddy', ds_eddy);

MeanFlow = struct('U_m', U_m, 'GradU_m', GradU_m, 'type', type);
end
