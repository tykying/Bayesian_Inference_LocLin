% Visualise concentration of particles
L = 1;

% Use saved statistics
Use_Saved_Stat = 0;

% TS Tracking Inference
TrajType_Save = 'ParticleConc';
%DataTypeString = '_FrequStat';
%nparticles = 2000;
Nx_cell_options = 16;
DGscheme = 'U0K0';
veloc_testcase = 'QGM2'

%load('/home/s1046972/opt/qgm2_particle/PART_TRAJ_nu300/traj_Npart10000_tIntv5yr_SampletIntv24hr.mat')
load('/home/s1046972/opt/qgm2_particle/PART_TRAJ/traj_Npart20000_tIntv8yr_SampletIntv24hr.mat')

% Read x,y,ts_list,L,ngrid
%%% Convert x, y from centimeters to meters
x = x/(100); y = y/(100);  L = L/(100);

%h = (ts_list(2)-ts_list(1))/(3600);
h = (ts_list(2)-ts_list(1));

time_scale = 3600*24;
space_scale = 1000;

SamplingInterval_factor_list = [1];

Nx_cell = 20;
Ny_cell = Nx_cell;
SpDis_ranges_min = [0, 0];
SpDis_ranges_max = [L, L]/space_scale;

[Mesh, RectMesh_Param] = setup_RectMeshStruct(Nx_cell, Ny_cell, SpDis_ranges_min, SpDis_ranges_max);
bin_size_x = Mesh(1).cell_dx(1);
bin_size_y = Mesh(1).cell_dx(2);

grid_vis = setup_grid_vis(Mesh, RectMesh_Param)

for SamplingInterval_factor_iterator = 1:length(SamplingInterval_factor_list)
    SamplingInterval_factor = SamplingInterval_factor_list(SamplingInterval_factor_iterator);
    SamplingParticle = 1;
    
    h = (ts_list(2) - ts_list(1));
    SamplingInterval = SamplingInterval_factor*h;
    
    TrajData_SubSampled = SubSample_Traj(x, y, ts_list, SamplingInterval, SamplingParticle);
    clear_vars = {'x', 'y', 'ts_list'};
    
    x_SubSampled = TrajData_SubSampled.x_SubSampled/space_scale;
    y_SubSampled = TrajData_SubSampled.y_SubSampled/space_scale;
    ts_list_SubSampled = TrajData_SubSampled.ts_list_SubSampled/time_scale;
    
    nparticles = size(x_SubSampled, 2);
    
%     diffx_SubSampled = diff(x_SubSampled, 1);
%     diffy_SubSampled = diff(y_SubSampled, 1);
%     
%     h_SubSampled = zeros(size(diffx_SubSampled));
%     for part = 1:size(h_SubSampled, 2)
%         h_SubSampled(:, part) = diff(ts_list_SubSampled);
%     end

    PartConc = zeros(length(ts_list_SubSampled), Nx_cell, Ny_cell);
    
    % Visualise Jumps distribution
    for t_ind = 1:length(ts_list_SubSampled)
        tn = ts_list_SubSampled(t_ind)-ts_list_SubSampled(1);
        
        % P.D.F. of Particle Concentrations
        i_CellInd = ceil(x_SubSampled(t_ind, :)/bin_size_x);
        j_CellInd = ceil(y_SubSampled(t_ind, :)/bin_size_y);
        k_CellInd = RectMesh_Param.cellij_to_k(i_CellInd, j_CellInd);
        
        [N_k, edges] = histcounts(k_CellInd, Nx_cell*Ny_cell);
        N_ij = convert_cellk_to_cellij(N_k', RectMesh_Param);
        PartConc_Inst = N_ij/nparticles;

        % Output for storage
        PartConc(t_ind, :,:) = PartConc_Inst;
    end
end


% Analysis of the time-averaged concentration field
SS_begin = round(length(ts_list_SubSampled)*0.6);  % use a different begin time
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

%% Visualisation
% Relative Derivation

f60 = figure(60)
set(f60, 'Position', [100, 100, 800, 800]);
subplot(2,1,1)
ContourLineValue = [-0.2:0.02:0.2];
contourf(grid_vis.X, grid_vis.Y, PartConc_TimeAverage_Deviation, ContourLineValue)
caxis([ContourLineValue(1), ContourLineValue(end)]);
xlim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
ylim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
colorbar;
pbaspect([1 1 1])
title('PartConc; All Data');


subplot(2,1,2)
ContourLineValue = [-0.2:0.02:0.2];
contourf(grid_vis.X, grid_vis.Y, PartConc_SSTimeAverage_Deviation, ContourLineValue)
caxis([ContourLineValue(1), ContourLineValue(end)]);
xlim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
ylim([RectMesh_Param.grid_x(1), RectMesh_Param.grid_x(end)]);
colorbar;
pbaspect([1 1 1])
title('Post steady-state');

close all;

% Make video of the P.D.F. of Jumps at different sampling interval
TrajType_Save = 'ParticleConc';
video = VideoWriter(['./Video_Generated/', TrajType_Save, '_', veloc_testcase, '_TimeEvol.avi'],'Uncompressed AVI');
video.FrameRate = 10
open(video)
vframe = figure(22)
colormap(jet(256))
set(vframe, 'Position', [100, 100, 800, 800]);
for SamplingInterval_factor_iterator = 1:length(SamplingInterval_factor_list)
    % Visualise Jumps distribution
    for t_ind = 1:length(ts_list_SubSampled)
        tn = ts_list_SubSampled(t_ind)-ts_list_SubSampled(1);
        
        % P.D.F. of Particle Concentrations
        %ContourLineValue = [0:0.0005:0.010];
        ContourLineValue = [-0.3:0.05:0.3];
        contourf(grid_vis.X, grid_vis.Y, squeeze(PartConc_RelDeviation(t_ind,:,:)), ContourLineValue)
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

close all;


% Visualisation for Mixing
% Linear concentration
TrajType_Save = 'ParticleMix';
video = VideoWriter(['./Video_Generated/', TrajType_Save, '_', veloc_testcase, '_TimeEvol.avi'],'Uncompressed AVI');
video.FrameRate = 10;
open(video)

f11 = figure(11);
set(f11, 'Position', [100, 100, 800, 800]);
colormap(jet(256))
c = linspace(1,25,nparticles);
c = fliplr(c);

for t_ind = 1:length(ts_list_SubSampled)
    tn = ts_list_SubSampled(t_ind) - ts_list_SubSampled(1);
    %plot_MeanFlowStreamline(space_scale);
    
    scatter(x_SubSampled(t_ind,:), y_SubSampled(t_ind,:), [], c, 'filled');
    
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
