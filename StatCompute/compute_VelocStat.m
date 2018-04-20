function [VelocStat] = compute_VelocStat(TrajData_Binned, VelocData, Mesh_Struct)
NCells = length(Mesh_Struct.Mesh);

%x = TrajData_Binned.x;
%y = TrajData_Binned.y;
%ts_list = TrajData_Binned.ts_list;
GridInd = TrajData_Binned.GridInd;
GridInd(end, :) = [];  % Ito: shift from t_n to t_n; Remove the last row for displacement

u = VelocData.u;
v = VelocData.v;

%assert(all(size(u) == (size(x) - [1, 0])));

% 4 Moments; 2 Component: u,v
Moments_eddies = zeros(NCells, 4, 2);

% %%% Only when spatial structure is involved
% x(end, :) = [];
% y(end, :) = [];
for cell_k = 1:NCells
    GridInd_inCell = (GridInd == cell_k);
    
%   %%% Only when spatial structure is involved
%    x_inCell = x(GridInd_inCell);
%    y_inCell = y(GridInd_inCell);
    
    u_inCell = u(GridInd_inCell);
    v_inCell = v(GridInd_inCell);
    
    U_inCell = [u_inCell, v_inCell];
    
    Moments_eddies(cell_k, 1, :) = mean(U_inCell, 1);
    Moments_eddies(cell_k, 2, :) = var(U_inCell, 1);
    Moments_eddies(cell_k, 3, :) = skewness(U_inCell, 1);
    Moments_eddies(cell_k, 4, :) = kurtosis(U_inCell, 1);
end

VelocStat = struct('Moments_eddies', Moments_eddies);