% Visualise concentration of particles
%% Loading And Initialise Data
% Load Trajectory Data
veloc_testcase = 'QGM2';
run('Script_LoadTrajData');

% Setup Rectangular Mesh
Nx_cell_options = 32;
run('Script_SetupRectGrid');

% Sub-Sampling Trajectory
SamplingInterval = h;
SamplingParticle = 1;
t_begin = ts_list(1) + 6 * (365*24*3600);
run('Script_SubSampling_Traj');

%% Calculate concentration of particles in rectangular bins
% Rescale the xample
x = x/space_scale;
y = y/space_scale;
ts_list = ts_list/time_scale;

[Nts, nparticles] = size(x);

%     diffx = diff(x, 1);
%     diffy = diff(y, 1);
%
%     h_SubS = zeros(size(diffx));
%     for part = 1:size(h_SubS, 2)
%         h_SubS(:, part) = diff(ts_list);
%     end

PartConc = zeros(Nts, Nx_cell, Ny_cell);

% Visualise Jumps distribution
for i_t = 1:Nts
    tn = ts_list(i_t)-ts_list(1);
    
    % P.D.F. of Particle Concentrations
    i_CellInd = ceil(x(i_t, :)/bin_sizes(1));
    j_CellInd = ceil(y(i_t, :)/bin_sizes(2));
    k_CellInd = RectMesh_Param.cellij_to_k(i_CellInd, j_CellInd);
    
    [N_k, edges] = histcounts(k_CellInd, Nx_cell*Ny_cell);
    N_ij = convert_cellk_to_cellij(N_k', RectMesh_Param);
    PartConc_Inst = N_ij/nparticles;
    
    % Output for storage
    PartConc(i_t, :,:) = PartConc_Inst;
end


% Analysis of the time-averaged concentration field
SS_begin = round(Nts*0.7);
PartConc_TimeAverage = squeeze(mean(PartConc, 1));
PartConc_SSTimeAverage = squeeze(mean(PartConc(SS_begin:end,:,:), 1));

PartConc_EvenSpread = (nparticles/(Nx_cell*Ny_cell))/nparticles;  % Assume equal grid size

PartConc_RelDeviation = (PartConc - PartConc_EvenSpread)/PartConc_EvenSpread;

PartConc_TimeAverage_Deviation = (PartConc_TimeAverage - PartConc_EvenSpread)/PartConc_EvenSpread;
PartConc_SSTimeAverage_Deviation = (PartConc_SSTimeAverage - PartConc_EvenSpread)/PartConc_EvenSpread;

PartConc_TimeAverage_SD = mean(abs(PartConc_TimeAverage_Deviation(:)))
PartConc_TimeAverage_SDmed = median(abs(PartConc_TimeAverage_Deviation(:)))

PartConc_SSTimeAverage_SD = mean(abs(PartConc_SSTimeAverage_Deviation(:)))
PartConc_SSTimeAverage_SDmed = median(abs(PartConc_SSTimeAverage_Deviation(:)))


% Refinement for easy visualisation
refinement_ratio = [2,2];

% Obtain the corresponding mesh
[Mesh_refined, RectMesh_Param_refined] = setup_RectMeshStruct(refinement_ratio(1)*Nx_cell, refinement_ratio(2)*Ny_cell, SpDis_ranges_min, SpDis_ranges_max);
grid_vis_refined = setup_grid_vis(Mesh_refined, RectMesh_Param_refined)


PartConc_RelDeviation_refined = zeros([size(PartConc_RelDeviation).*[1, refinement_ratio]]);
for i_t = 1:size(PartConc_RelDeviation_refined, 1)
    PartConc_RelDeviation_refined(i_t,:,:) = griddata_interpolation(squeeze(PartConc_RelDeviation(i_t,:,:)), refinement_ratio);
end

PartConc_TimeAverage_Deviation_refined = griddata_interpolation(PartConc_TimeAverage_Deviation, refinement_ratio);
PartConc_SSTimeAverage_Deviation_refined = griddata_interpolation(PartConc_SSTimeAverage_Deviation, refinement_ratio);


stop
%% Visualisation
npart_scale_str = 1000;


% Relative Derivation
f60 = figure(60)
set(f60, 'Position', [100, 100, 800, 800]);
colormap(jet(256))
subplot(2,1,1)
ContourLineValue = [-0.10:0.02:0.10];
% contourf(grid_vis.X, grid_vis.Y, PartConc_TimeAverage_Deviation, ContourLineValue)
contourf(grid_vis_refined.X, grid_vis_refined.Y, PartConc_TimeAverage_Deviation_refined, ContourLineValue)
caxis([ContourLineValue(1), ContourLineValue(end)]);
xlim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
ylim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
colorbar;
pbaspect([1 1 1])
title('PartConc; All Data');


subplot(2,1,2)
ContourLineValue = [-0.10:0.02:0.10];
% contourf(grid_vis.X, grid_vis.Y, PartConc_SSTimeAverage_Deviation, ContourLineValue)
contourf(grid_vis_refined.X, grid_vis_refined.Y, PartConc_SSTimeAverage_Deviation_refined, ContourLineValue)
caxis([ContourLineValue(1), ContourLineValue(end)]);
xlim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
ylim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
colorbar;
pbaspect([1 1 1])
title('Post steady-state');

close all;

% Make video of the P.D.F. of Jumps at different sampling interval
TrajType_Save = 'ParticleConc';
video = VideoWriter(['./Video_Generated/', TrajType_Save, '_', veloc_testcase, '_npart', num2str(nparticles/npart_scale_str), 'k','_Idx',num2str(Nx_cell_options), '_TimeEvol.avi'],'Uncompressed AVI');
video.FrameRate = 10
open(video)
vframe = figure(22)
colormap(jet(256))
set(vframe, 'Position', [100, 100, 800, 800]);
for SamplingInterval_factor_iterator = 1:length(SamplingInterval_factor_list)   
    % Visualise Jumps distribution
    for i_t = 1:Nts
        tn = ts_list(i_t)-ts_list(1);
        
        % P.D.F. of Particle Concentrations
        %ContourLineValue = [0:0.0005:0.010];
        ContourLineValue = [-0.3:0.05:0.3];
        %contourf(grid_vis.X, grid_vis.Y, squeeze(PartConc_RelDeviation(i_t,:,:)), ContourLineValue)
       	contourf(grid_vis_refined.X, grid_vis_refined.Y, squeeze(PartConc_RelDeviation_refined(i_t,:,:)), ContourLineValue)
        caxis([ContourLineValue(1), ContourLineValue(end)]);
        xlim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
        ylim([RectMesh_Param.grid_y(1), RectMesh_Param.grid_y(end)]);

        title(['Particle Concentration on day ', num2str(tn)]);
        colorbar;
        pbaspect([1 1 1])
        
        drawnow();
        frame = getframe(vframe);
        writeVideo(video,frame);
        
        clf;      
    end
end
close(video)

PartConc_RelDeviation_refined = griddata_interpolation(squeeze(PartConc_RelDeviation(i_t,:,:)), [2,2]);

close all;

stop

% Visualisation for Mixing
% Linear concentration
TrajType_Save = 'ParticleMix';
video = VideoWriter(['./Video_Generated/', TrajType_Save, '_', veloc_testcase, '_', num2str(nparticles/npart_scale_str), 'k_TimeEvol.avi'],'Uncompressed AVI');
video.FrameRate = 10;
open(video)

f11 = figure(11);
set(f11, 'Position', [100, 100, 800, 800]);
colormap(jet(256))
c = linspace(1,25,nparticles);
c = fliplr(c);

for i_t = 1:Nts
    tn = ts_list(i_t) - ts_list(1);
    %plot_MeanFlowStreamline(space_scale);
    
    scatter(x(i_t,:), y(i_t,:), [], c, 'filled');
    
    tn_string = sprintf('%4.2f', tn);
    
    xlim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
    ylim([RectMesh_Param.grid_y(1), RectMesh_Param.grid_y(end)]);
    %title(['time = ', tn_string, ' months'])
    %xlabel('x (in km)');
    %xlabel('y (in km)');
    title(['day = ', tn_string])
    xlabel('x (in km)');
    ylabel('y (in km)');
    pbaspect([1 1 1])
    alpha(0.65)
    
    drawnow();
    frame = getframe(f11);
    writeVideo(video,frame);
    
    clf
end
close(video)
