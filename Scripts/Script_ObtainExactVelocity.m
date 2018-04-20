addpath('../PlotUtilities/')
addpath('../Utilities/')
addpath('../StructConversion/')

% Assume Exact Field came from inference
% This script extract the inferred velocity field
Set = 10;
if Set == 10
    Layer = 2;
    Nx_cell_ARG = 32;
    SamplingInterval_vis = 1;
    DS_rate = 1;
    Nsteps_pc = 4096;
    nTrial_List = [1:1];
    
    veloc_Profile = 'QGM2_DSpart';
    InfScheme = 'Local';
elseif Set == 11
    Layer = 2;
    Nx_cell_ARG = 16;
    SamplingInterval_vis = 1;
    DS_rate = 1;
    Nsteps_pc = 4096;
    nTrial_List = [1:1];
    
    veloc_Profile = 'QGM2_nReal_DSTemp';
    InfScheme = 'Local';
else
    Layer = 2;
    Nx_cell_ARG = 64;
    SamplingInterval_vis = 1;
    DS_rate = 4;
    Nsteps_pc = 4096*2;
    nTrial_List = [1:1];
    
    veloc_Profile = 'QGM2_nReal_DSTemp';
    InfScheme = 'Local';
end

for nTrial = nTrial_List
DataTypeString = 'BinnedTrajData';
run('../Scripts/Script_Filenames');
BinnedTrajData = load(filename_BTD);

DataTypeString = 'BayesSampleStat';
run('../Scripts/Script_Filenames');
BayesSampleStat = load(filename_BSS);

Mesh_Struct = BayesSampleStat.Mesh_Struct;
RectMesh_Param = Mesh_Struct.RectMesh_Param;
Mesh = Mesh_Struct.Mesh;


% PWC Velocity Fields
diffx_moments = BinnedTrajData.TrajJumps_DA.TrajJumps_MomentGlobal.diffx;
diffy_moments = BinnedTrajData.TrajJumps_DA.TrajJumps_MomentGlobal.diffy;

% Assumed uniform time-interval
h_moments = BinnedTrajData.TrajJumps_DA.TrajJumps_MomentGlobal.h;

u_pwc = diffx_moments(:,1)./h_moments(:,1);
v_pwc = diffy_moments(:,1)./h_moments(:,1);

grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
u_pwc_ij = convert_cellk_to_cellij( u_pwc, RectMesh_Param );
v_pwc_ij = convert_cellk_to_cellij( v_pwc, RectMesh_Param );

figure(Nx_cell_ARG+1+Nsteps_pc)
contourf(grid_vis.X, grid_vis.Y, u_pwc_ij)
caxis([-0.3, 0.3])
%plot_MeanFlowStreamline(gca, grid_vis, u_pwc_ij, v_pwc_ij)

figure(Nx_cell_ARG+2+Nsteps_pc)
contourf(grid_vis.X, grid_vis.Y, v_pwc_ij)
caxis([-0.3, 0.3])
%plot_MeanFlowStreamline(gca, grid_vis, u_pwc_ij, v_pwc_ij)


% PWL Velocity Fields
MCMC_NSample = BayesSampleStat.MCMC_Stat.MCMC_NSample;
theta_mean = BayesSampleStat.theta_Stat.theta_sum/MCMC_NSample;
thetasq_mean = BayesSampleStat.theta_Stat.thetasq_sum/MCMC_NSample;


[MeanFields_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_mean, Mesh_Struct);
u_pwl_ij = MeanFields_f.u;
v_pwl_ij = MeanFields_f.v;

figure(Nx_cell_ARG+10+Nsteps_pc)
contourf(grid_vis_f.X, grid_vis_f.Y, u_pwl_ij)
caxis([-0.3, 0.3])
%plot_MeanFlowStreamline(gca, grid_vis_f, u_pwl_ij, v_pwl_ij)

figure(Nx_cell_ARG+20+Nsteps_pc)
contourf(grid_vis_f.X, grid_vis_f.Y, v_pwl_ij)
caxis([-0.3, 0.3])
%plot_MeanFlowStreamline(gca, grid_vis_f, u_pwl_ij, v_pwl_ij)


% Smoothen u_pwl_ij
u_pwl_ij_sm = u_pwl_ij;
v_pwl_ij_sm = v_pwl_ij;
for i=2:(size(u_pwl_ij, 1)-1)
    for j=2:(size(u_pwl_ij, 2)-1)
        u_pwl_ij_sm(i, j) = (sum(u_pwl_ij(i-1:i+1,j))+sum(u_pwl_ij(i,j-1:j+1)))/6;
        v_pwl_ij_sm(i, j) = (sum(v_pwl_ij(i-1:i+1,j))+sum(v_pwl_ij(i,j-1:j+1)))/6;
    end
end

figure(Nx_cell_ARG+100+Nsteps_pc)
contourf(grid_vis_f.X, grid_vis_f.Y, u_pwl_ij_sm)
caxis([-0.3, 0.3])
plot_MeanFlowStreamline(gca, grid_vis_f, u_pwl_ij, v_pwl_ij);

figure(Nx_cell_ARG+200+Nsteps_pc)
contourf(grid_vis_f.X, grid_vis_f.Y, v_pwl_ij_sm)
caxis([-0.3, 0.3])
plot_MeanFlowStreamline(gca, grid_vis_f, u_pwl_ij, v_pwl_ij);

end