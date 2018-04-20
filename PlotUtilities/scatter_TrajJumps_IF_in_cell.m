% Show distribution of jumps
%% Usual set of parameters
% Parameters in sliders (Integer-valued, from 1 to max)
Layer_List = [2];
Nx_cell_ARG_List = [32  16];
SamplingInterval_vis_List = [8, 32, 16, 1, 2, 4, 64];

% Configurations for video
Nlevel = 1;
veloc_Profile = 'QGM2_nReal_DSTemp';
nReal = 60;
DS_rate = 4;   % DownSampler rate
traj_fullpath_old = '';

for Layer = Layer_List
for SamplingInterval_vis = SamplingInterval_vis_List
for Nx_cell_ARG = Nx_cell_ARG_List
    assert(nReal == 60);
    assert(DS_rate == 4);

    % Load Data
    run('./Scripts/Script_LoadTrajData');
    disp('Finished Loading Raw Trajectory Data.')
    
    % Sub-sampling Trajectories
    run('./Scripts/Script_SubSampling_Traj');
    disp('Finished SubSampling Trajectories.')
    
    % Set-up Mesh
    run('./Scripts/Script_SetupRectMesh');
    disp('Finished Setting Up RectMesh.')
    
    % Sort TrajJumps
    [GridInd, alphax, alphay] = positioning_RectMesh(TrajData.x, TrajData.y, RectMesh_Param);
    TrajPosition = struct('GridInd', GridInd, ...
        'alphax', alphax, 'alphay', alphay );
    disp('Finished positioning and computing alpha.')
    
    [TrajJumps, NJumps_cell] = sort_TrajJumps_Ito_MultipleSteps(TrajData, TrajPosition, Mesh, Nlevel);
    Dt_sortJumps= toc;
    assert(all(NJumps_cell >= 5))
    disp('Finished Sorting TrajJumps.')
    
    % Video Production
    DataTypeString = 'ParticleBinnedEvol';
    run('./Scripts/Script_Filenames');
    [status, msg, msgID] = mkdir(filepath);
    
    %video = VideoWriter(filename_PBE, 'Uncompressed AVI');
    video = VideoWriter(filename_PBE, 'Motion JPEG AVI');

    % Below: Only for QGM2_nReal_DSTemp
    assert(strcmp(veloc_Profile, 'QGM2_nReal_DSTemp'));
    
    % ad-hoc for QGM2_nReal_DSTemp
    nparticles_pc = length(TrajJumps(1).x0);
    nparticles_pcpr = nparticles_pc/nReal;    
    
    video.FrameRate = 1;
    open(video)
    vframe = figure(22)
    colormap(jet(256))
    set(vframe, 'Position', [100, 100, 1280, 720]);
    
    for cell_inv_IND = 1:length(TrajJumps)
        cell_inv_X = Mesh(cell_inv_IND).cell_centre;
        %     for part_Ind_inv = 1:( 2^(log2(nparticles_pcpr)/2) + 2^log2(nparticles_pcpr)/8 ):nparticles_pcpr
        %         part_Ind_List = part_Ind_inv:nparticles_pcpr:nparticles_pc;
        part_Ind_List = 1:nparticles_pc;
        
        clf
        TrajJumps_inv = TrajJumps(cell_inv_IND);
        part_x_I = (TrajJumps_inv.x0 + cell_inv_X(1))/space_scale;
        part_y_I = (TrajJumps_inv.y0 + cell_inv_X(2))/space_scale;
        
        part_x_F = (TrajJumps_inv.x0 + TrajJumps_inv.diffx + cell_inv_X(1))/space_scale;
        part_y_F = (TrajJumps_inv.y0 + TrajJumps_inv.diffy + cell_inv_X(2))/space_scale;
        
        scatter(part_x_F(part_Ind_List), part_y_F(part_Ind_List), 'r.');
        hold on
        scatter(cell_inv_X(1)/space_scale, cell_inv_X(2)/space_scale, 'b');
%         scatter(part_x_I(part_Ind_List), part_y_I(part_Ind_List), 'b.');
        
        x_range = [RectMesh_Param.SpDis_ranges_min(1), RectMesh_Param.SpDis_ranges_max(1)]/space_scale;
        y_range = [RectMesh_Param.SpDis_ranges_min(2), RectMesh_Param.SpDis_ranges_max(2)]/space_scale;
        
        xlim(x_range);
        ylim(y_range);
        axis square
        
        cell_inv_IND_ij = RectMesh_Param.cellk_to_ij(cell_inv_IND);
        title(sprintf('In cell (%d, %d); h= %2.1f days; NPart= %d; nReal=%d', cell_inv_IND_ij(1), cell_inv_IND_ij(2), h/time_scale, nparticles_pc, nReal))

        set(gca,'xtick',linspace(x_range(1), x_range(2), Nx_cell_ARG+1));
        set(gca,'ytick',linspace(y_range(1), y_range(2), Nx_cell_ARG+1));
        
        xlabel('x (in km)');
        ylabel('y (in km)');

        grid on

        drawnow()
        frame = getframe(vframe);
        writeVideo(video,frame);
        %     end
    end
    close(video)

end
end
end