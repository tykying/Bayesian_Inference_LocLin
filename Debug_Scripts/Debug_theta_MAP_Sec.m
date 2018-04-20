%% Visualise theta_MAP_Sec after the long initial MCMC Chain
addpath('./PlotUtilities/')
grid_vis = setup_grid_vis(Mesh, RectMesh_Param);

% Per sessions
for Ses = 1:16:NSessions
    theta_MAP_Sec = theta_MAP_LIST(:, :, Ses);
    
    
    figure(200+Ses)
    theta_MAP_Sec_ij = convert_cellk_to_cellij( theta_MAP_Sec, RectMesh_Param );
    ax_sp = plot_dataij_k(grid_vis, theta_MAP_Sec_ij);
end

% Global MAP

% Manuelly remove a cell
cell_k_RM = RectMesh_Param.cellij_to_k(6, 10);

theta_MAP_Global(cell_k_RM, :) = NaN;


figure(300+0)
theta_MAP_Global_ij = convert_cellk_to_cellij( theta_MAP_Global, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_MAP_Global_ij);


%% For Picked Theta
figure(400)
theta_Init_Picked_ij = convert_cellk_to_cellij( theta_Init_Picked, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_Init_Picked_ij);

figure(500)
theta_SSD_Picked_ij = convert_cellk_to_cellij( theta_SSD_Picked, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_SSD_Picked_ij);

figure(500)
theta_SSD_Picked_ij = convert_cellk_to_cellij( theta_SSD_Picked, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_SSD_Picked_ij);
