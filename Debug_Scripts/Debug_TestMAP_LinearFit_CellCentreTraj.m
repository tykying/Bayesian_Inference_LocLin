% Script to test the
%% Obtain the MAP theta
logPost_ts = theta_Stat.logPost_ts;

logPost_MAX = max(logPost_ts, [], 2);
theta_MAP = zeros(length(Mesh), NVars);
MAP_SampleInd = zeros(length(Mesh), 1);
for cell_k = 1:length(Mesh)
    MAP_SampleInd(cell_k) = find(logPost_ts(cell_k, :) == logPost_MAX(cell_k), 1);
    
    theta_MAP(cell_k, :) = squeeze( theta_Stat.theta_store(cell_k,:,MAP_SampleInd(cell_k)) );
end

theta_MAP = convert_theta_U1K0_to_theta_KPolar(theta_MAP);

theta_LinFit = linear_fit_TrajJumps(TrajJumps);
IncompressiblityCheck = sum(theta_LinFit(:, [3, 6]), 2)./max(abs(theta_LinFit(:, [3, 6])), [], 2);



%% FOR DEBUGGING
figure(10)
logPost_ts_vis = logPost_ts(1:2:end, 1:8:end);
plot((logPost_ts_vis'));

addpath('./PlotUtilities/')
grid_vis = setup_grid_vis(Mesh, RectMesh_Param);

figure(1)
theta_MAP_ij = convert_cellk_to_cellij( theta_MAP, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_MAP_ij);

figure(2)
AccRatio = MCMC_Stat.acceptance_ratio;
AccRatio_ij = convert_cellk_to_cellij( AccRatio, RectMesh_Param );

ax_sp = plot_dataij_k(grid_vis, AccRatio_ij);

figure(3)
theta_MAP_FIL = theta_MAP;
theta_MAP_FIL(AccRatio > 0.4) = NaN;

theta_MAP_FIL_ij = convert_cellk_to_cellij( theta_MAP_FIL, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_MAP_FIL_ij);

figure(4)
MAP_SampleInd_ij = convert_cellk_to_cellij( MAP_SampleInd, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, MAP_SampleInd_ij);

figure(5)
theta_Scale_ij = convert_cellk_to_cellij( theta_Stat.theta_Scale, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_Scale_ij);

figure(6)
theta_LinFit_ij = convert_cellk_to_cellij( theta_LinFit, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_LinFit_ij);

figure(7)
IncompressiblityCheck_ij = convert_cellk_to_cellij( IncompressiblityCheck, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, IncompressiblityCheck_ij);


% For given cells, using the MAP field and compare with the particle
% displacement
cellij_to_k = RectMesh_Param.cellij_to_k;
cellk_List = [cellij_to_k(4, 11), cellij_to_k(6, 9), cellij_to_k(1, 10), cellij_to_k(11, 6) ...
    cellij_to_k(3, 3), cellij_to_k(9, 9)];
%%
cellk_List = [cellij_to_k(6, 3), cellij_to_k(8, 8), cellij_to_k(7, 16)];
for cell_k = cellk_List
    cellij = RectMesh_Param.cellk_to_ij(cell_k);
    theta_cell_k = theta_MAP(cell_k, :);
    
    figure(cell_k + 10000);
    [ax_quiver] = visualise_LinearVelocityFld(theta_cell_k, bin_sizes, SamplingInterval_vis)
    title(sprintf('cell (%d, %d)', cellij(1), cellij(2)) );
    hold on
    x0 = TrajJumps(cell_k).x0;
    y0 = TrajJumps(cell_k).y0;
    xf = TrajJumps(cell_k).x0 + TrajJumps(cell_k).diffx;
    yf = TrajJumps(cell_k).y0 + TrajJumps(cell_k).diffy;
    %scatter(x0, y0, 'r')
    %scatter(xf, yf, 'k')
end


% Round up order
round(log10(1.0001e-4)+0.5)