addpath('../PlotUtilities/')
addpath('../Scripts/')
addpath('./PlotUtilities/')
addpath('./Scripts/')



adhoc = 1;
if contains(RunProfile, 'QGM2')
    %% Cells of the interest
    if Nx_cell_ARG == 16
        %cell_ij_list = { [1, 9], [3, 9], [5, 10], [10, 10]; ...
        %    [14, 3], [15, 14], [3, 10], [3, 8]};
        
        %ij_list = [2, 6, 10, 14];
        ij_list = [3, 7, 12];
        ij_list = [3, 9, 13];

        x_list = ij_list;
        y_list = flip(ij_list);
        cell_ij_list = cell(length(ij_list), length(ij_list));

        % Arrange in the manner same as plot output
        for k_i = 1:length(x_list)
            for k_j = 1:length(y_list)
                cell_ij_list{k_j, k_i} = [x_list(k_i), y_list(k_j)];
            end
        end
        
        
        % For presentation: pick a cell to focus on
        cell_ij_Success = [11, 7];
        cell_ij_Fail = [8, 11];
        
        cell_ij_vis = [cell_ij_Success; cell_ij
            ail];
        
        cellij_to_k = RectMesh_Param.cellij_to_k;
        cellk_to_ij = RectMesh_Param.cellk_to_ij;
        
        cell_k_vis = zeros(1, size(cell_ij_vis, 1));
        for i = 1:size(cell_ij_vis, 1)
            cell_k_vis(i) = cellij_to_k(cell_ij_vis(i, 1), cell_ij_vis(i, 2));
        end
    end
    
    
end    

if (contains(RunProfile, 'QGM2') && contains(RunProfile, 'Plot_Traj'))
    %% 
    NParts_vis = 64;
    rng(NParts_vis);
    L_conv = 1000;

    NParts_List = sort(round(rand(1, NParts_vis)*size(x, 2)));
    x_vis = x(:, NParts_List)/L_conv;
    y_vis = y(:, NParts_List)/L_conv;
    
    FigObj = figure('Name', 'Fields_TrajVis', 'NumberTitle','on', 'pos', [200 200 800 800]);
    ax_cur_cur = gca;
    img_str = plot(ax_cur_cur, x_vis, y_vis);
    
    grid(ax_cur_cur, 'on')
    Grid_min = RectMesh_Param.SpDis_ranges_min/L_conv;
    Grid_max_cur = RectMesh_Param.SpDis_ranges_max_cur/L_conv;
    xlim(ax_cur_cur, [Grid_min(1), Grid_max_cur(1)]);
    ylim(ax_cur_cur, [Grid_min(1), Grid_max_cur(1)]);
    
    ax_cur_cur.XTick = linspace(Grid_min(1), Grid_max_cur(1), RectMesh_Param.Nx_cell+1);
    ax_cur_cur.YTick = linspace(Grid_min(2), Grid_max_cur(2), RectMesh_Param.Ny_cell+1);

    for i = 1:length(ax_cur_cur.XTick)
        if mod(i, 4) ~= 0
            ax_cur_cur.XTickLabels{i} = '';
        end
        if mod(i-1, 2) ~= 0
            ax_cur_cur.YTickLabels{i} = ''; 
        end
    end
    xlabel(ax_cur_cur, '$x$ [km]');
    ylabel(ax_cur_cur, '$y$ [km]');
    
    pbaspect(ax_cur_cur, [1 1 1])
    %colorbar;
    run('Script_ax_curesConfig.m')
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_NParts', num2str(NParts_vis)];
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_TrajVis', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
    % Draw rectangles to highlight a cell Additionally
    cell_ij_List = [cell_ij_Success; cell_ij_Fail];
    for i = 1:size(cell_ij_List, 1)
        Rect_ij = cell_ij_List(i, :);
        
        grid_x = RectMesh_Param.grid_x;
        grid_y = RectMesh_Param.grid_y;
        pos_ll = [grid_x(Rect_ij(1)), grid_y(Rect_ij(2))]/L_conv;
        pos_dxy = RectMesh_Param.bin_sizes/L_conv;
        Rect_pos = [pos_ll, pos_dxy];
        hold(ax_cur_cur, 'on')
        rectangle(ax_cur_cur, 'Position', Rect_pos, 'LineWidth',5)
    end
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_NParts', num2str(NParts_vis)];
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_TrajVisWithRect', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end

if (contains(RunProfile, 'QGM2') && contains(RunProfile, 'Histograms'))
    %% Pick the cells to focus on
    cell_ij_vis = [cell_ij_Success; cell_ij_Fail];
    
    cellij_to_k = RectMesh_Param.cellij_to_k;
    cellk_to_ij = RectMesh_Param.cellk_to_ij;
    
    cell_k_vis = zeros(1, size(cell_ij_vis, 1));
    for i = 1:size(cell_ij_vis, 1)
        cell_k_vis(i) = cellij_to_k(cell_ij_vis(i, 1), cell_ij_vis(i, 2));
    end
    
    nTrial_str = num2str(nTrial);
    
    MCMC_NSample = MCMC_Stat.MCMC_NSample;
    theta_store = theta_Stat.theta_store;
    
    [Ncells, NVars, MCMC_NSample] = size(theta_store);

    tic;
    NSessions = 1;
    if NSessions > 1
        theta_Fields_Ses = cell(1,NSessions);
        for Ses = 1:NSessions
            Ses_bgn = round(MCMC_NSample*(Ses-1)/NSessions)+1;
            Ses_end = round(MCMC_NSample*Ses/NSessions);
            fprintf('Ses: %d, in samples (%d, %d). \n', Ses, Ses_bgn, Ses_end);
            
            theta_Fields_Ses{Ses} = convert_theta_U1K0_to_theta_KPolar(theta_store(cell_k_vis, :, Ses_bgn:Ses_end));
            
            if Ses > 1
                theta_Fields = cat(3, theta_Fields, theta_Fields_Ses{Ses});
            else
                theta_Fields = theta_Fields_Ses{Ses};
            end
        end
    else
        theta_Fields = convert_theta_U1K0_to_theta_KPolar(theta_store(cell_k_vis, :, :));
    end
    DT_conKpolar = toc;
    fprintf('Finsihed Sesing MC in %f sec', DT_conKpolar);
    
    [K11, K22, K12] = Kpolar_to_Kcart_vectorised(theta_Fields(:,7,:), theta_Fields(:,8,:), theta_Fields(:,9,:));
    theta_Fields(:,7:9,:) = [K11, K22, K12];
    
    % Visualisation
    NCells_vis = length(cell_k_vis);
    
    %% U,V: on same plot
    
    close all;
    FigObj_List = {};
    FigObj_List{1} = figure('Name', FigObj_Name ,'NumberTitle','on', 'pos', [200 200 800 600]);
    FigObj = FigObj_List{1};
    for cell_ind = 1:NCells_vis
        cell_k = cell_k_vis(cell_ind);
        cellij = cellk_to_ij(cell_k_vis(cell_ind));
        
        FigObj_Name = sprintf('i: %d, j: %d; k: %d', cellij(1), cellij(2), cell_k)
        %FigObj_List{cell_ind} = figure('Name', FigObj_Name ,'NumberTitle','on', 'pos', [200 200 800 800]);
        
        %FigObj = FigObj_List{cell_ind};
        histfig = cell(NCells_vis, NVars);

        cell_k = cellij_to_k(cellij(1), cellij(2));
       
        
        Data_Stored = {theta_Fields(cell_ind, 1, :), theta_Fields(cell_ind, 2, :), {}; ...
                       theta_Fields(cell_ind, 4, :), theta_Fields(cell_ind, 5, :), theta_Fields(cell_ind, 6, :); ...
                       theta_Fields(cell_ind, 7, :), theta_Fields(cell_ind, 8, :), theta_Fields(cell_ind, 9, :); };
                   
                   
        title_label = {'$b_1$ [m/s]', '$b_2$ [m/s]', ''; ...
                       '$A_{11}$ [s$^{-1}$]', '$A_{12}$ [s$^{-1}$]', '$A_{21}$ [s$^{-1}$]',; ...
                       '$\kappa_{11} [m^2/s]$', '$\kappa_{22} [m^2/s]$', '$\kappa_{12} [m^2/s]$';};
                   
        title_label = {'$u$ [m/s]', '$v$ [m/s]', ''; ...
                       '$\partial_x u$ [s$^{-1}$]', '$\partial_y u$ [s$^{-1}$]', '$\partial_x v$ [s$^{-1}$]',; ...
                       '$\kappa_{11} [m^2/s]$', '$\kappa_{22} [m^2/s]$', '$\kappa_{12} [m^2/s]$';};    
                   
                   
        Data_Stored = {theta_Fields(cell_ind, 1, :), theta_Fields(cell_ind, 2, :), {}; ...
                       theta_Fields(cell_ind, 3, :), theta_Fields(cell_ind, 4, :), theta_Fields(cell_ind, 5, :); };

        xlim_Stored = {[-0.0125, 0.020], [-0.0135, 0.0120], {}; ...
                       {},  {}, {} };
                   
        title_label = {'$u$ [m/s]', '$v$ [m/s]', ''; ...
                       '$\partial_x u$ [s$^{-1}$]', '$\partial_y u$ [s$^{-1}$]', '$\partial_x v$ [s$^{-1}$]';};
                   
        % Transpose to get the right ordering
        Data_Stored = Data_Stored';
        title_label = title_label';
        xlim_Stored = xlim_Stored';
        
        ax_suplots = cell(size(Data_Stored, 1), size(Data_Stored, 2));
        for cmp = 1:length(Data_Stored(:))
            if isempty(Data_Stored{cmp}) ~= 1
                ax_suplots{cmp} = subplot(size(Data_Stored', 1), size(Data_Stored', 2), cmp);
                ax_cur = gca;
                hold(ax_cur, 'on');
                
                Data_Plot = squeeze(Data_Stored{cmp});
                histfig{cell_ind, cmp} = histogram(ax_cur, Data_Plot, 24);
                histfig{cell_ind, cmp}.Normalization = 'pdf';
                title(ax_cur, title_label{cmp});
                if isempty(xlim_Stored{cmp}) ~= 1
                    xlim(xlim_Stored{cmp});
                end
                
                %pbaspect(ax_cur, [1.1 1 1])
                set(ax_cur.Title, 'FontSize', 22, 'Interpreter', 'latex')
                set(ax_cur.XAxis, 'FontSize', 13.5)
                %set(ax_cur,'ytick',[])

                %run('Script_ax_curesConfig.m')
            end
        end
        
        for subplot_id = [1, 4] % first column
            ylabel(ax_suplots{subplot_id}, 'p.d.f.')
            set(ax_suplots{subplot_id}.YLabel, 'FontSize', 20, 'Interpreter', 'latex')
        end
        
        %suptitle(sprintf('cell (%d, %d)', cellij(1), cellij(2)))
    end
    legend({'Outside jet', 'On jet'}, 'position', [0.70 0.70 0.2 0.1], 'FontSize', 20, 'Interpreter', 'latex');
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_slides', '_h', num2str(SamplingInterval_vis)];
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_HistUV', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
    %% Diffusivity: on different subplots
    close all;
    FigObj_List = {};
    FigObj_List{1} = figure('Name', FigObj_Name ,'NumberTitle','on', 'pos', [200 200 800 600]);
    FigObj = FigObj_List{1};
    
    cell_k = cell_k_vis(cell_ind);
    cellij = cellk_to_ij(cell_k_vis(cell_ind));
    
    FigObj_Name = sprintf('i: %d, j: %d; k: %d', cellij(1), cellij(2), cell_k)
    %FigObj_List{cell_ind} = figure('Name', FigObj_Name ,'NumberTitle','on', 'pos', [200 200 800 800]);
    
    %FigObj = FigObj_List{cell_ind};
    histfig = cell(NCells_vis, NVars);
       
    
    Data_Stored = {theta_Fields(1, 7, :), theta_Fields(1, 8, :), theta_Fields(1, 9, :); 
                    theta_Fields(2, 7, :), theta_Fields(2, 8, :), theta_Fields(2, 9, :); };
       
    title_label = {'$\kappa_{11} [m^2/s]$', '$\kappa_{22} [m^2/s]$', '$\kappa_{12} [m^2/s]$';
        '$\kappa_{11} [m^2/s]$', '$\kappa_{22} [m^2/s]$', '$\kappa_{12} [m^2/s]$';};
    
    xlim_Stored = {[0, 250], [0, 350], [-125, 125]; [0, 17500], [0, 7500], [-5000, 5000] };
    
    % Transpose to get the right ordering
    Data_Stored = Data_Stored';
    title_label = title_label';
    xlim_Stored = xlim_Stored';
    
    ax_suplots = cell(size(Data_Stored, 1), size(Data_Stored, 2));
    for cmp = 1:length(Data_Stored(:))
        if isempty(Data_Stored{cmp}) ~= 1
            ax_suplots{cmp} = subplot(size(Data_Stored', 1), size(Data_Stored', 2), cmp);
            ax_cur = gca;
            hold(ax_cur, 'on');
            
            Data_Plot = squeeze(Data_Stored{cmp});
            histfig{cell_ind, cmp} = histogram(ax_cur, Data_Plot, 24);
            histfig{cell_ind, cmp}.Normalization = 'pdf';
            title(ax_cur, title_label{cmp});
            if isempty(xlim_Stored{cmp}) ~= 1
                xlim(xlim_Stored{cmp});
            end
            
            %pbaspect(ax_cur, [1.1 1 1])
            set(ax_cur.Title, 'FontSize', 22, 'Interpreter', 'latex')
            set(ax_cur.XAxis, 'FontSize', 13.5)
            %set(ax_cur,'ytick',[])
            
            %run('Script_ax_curesConfig.m')
        end
    end
    
    for subplot_id = [1] % first column
        ylabel(ax_suplots{subplot_id}, 'p.d.f. [Outside Jet]')
        set(ax_suplots{subplot_id}.YLabel, 'FontSize', 20, 'Interpreter', 'latex')
    end
    
    for subplot_id = [4] % first column
        ylabel(ax_suplots{subplot_id}, 'p.d.f. [On Jet]')
        set(ax_suplots{subplot_id}.YLabel, 'FontSize', 20, 'Interpreter', 'latex')
    end
    
    for cmp = [4:6] % first column
        histfig{cell_ind, cmp}.FaceColor = [0.8500    0.3250    0.0980];
    end
    
    %suptitle(sprintf('cell (%d, %d)', cellij(1), cellij(2)))
    %legend({'Outside jet', 'On jet'}, 'position', [0.70 0.70 0.2 0.1], 'FontSize', 20, 'Interpreter', 'latex');
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_slides', '_h', num2str(SamplingInterval_vis)];
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_HistK', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end

if (contains(RunProfile, 'QGM2') && contains(RunProfile, 'UVK_contourf'))
    %% Pick the cells to focus on
    nTrial_str = num2str(nTrial);
    
    MCMC_Stat = MCMC_Stat_List{1};
    theta_Stat = theta_Stat_List{1};
    
    MCMC_NSample = MCMC_Stat.MCMC_NSample;
    theta_store = theta_Stat.theta_store;
    
    [Ncells, NVars, MCMC_NSample] = size(theta_store);

    % Obtain MAP
    logPost_ts = theta_Stat.logPost_ts;
    for cell_k = 1:Ncells
        if mod(log2(cell_k), 4) == 0
            fprintf('Working on %d. \n', cell_k);
        end
        logPost_ts_cell_k = logPost_ts(cell_k, :);
        
        SampInd = find(logPost_ts_cell_k == max(logPost_ts_cell_k(:)), 1);
        
        theta_MAP_cell_k = squeeze(theta_store(cell_k, :, SampInd));
        
        theta_MAP(cell_k, :) = theta_MAP_cell_k;
    end
    
    theta_MAP = convert_theta_U1K0_to_theta_KPolar(theta_MAP);
    [K11, K22, K12] = Kpolar_to_Kcart_vectorised(theta_MAP(:,7,:), theta_MAP(:,8,:), theta_MAP(:,9,:));
    theta_MAP(:,7:9,:) = [K11, K22, K12];
    
        
    [theta_MAP_ij_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_MAP, Mesh_Struct);
    
    
    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);

    theta_MAP_ij = convert_cellk_to_cellij(theta_MAP, RectMesh_Param);
%%

    SamplingInterval_vis = 1;
    Fig_outputpath_lbl = {'_U', '_V'};
    for lbl = 1:2
        FigObj = figure('Name', 'Fields_UV', 'NumberTitle','on', 'pos', [200 200 800 800]);
        FigObj_Name = sprintf('U: h= %d days', SamplingInterval_vis)
        ax_cur = gca;
        
        % Happen that 1,2 are also velocity component
        dataij = squeeze(theta_MAP_ij(:,:,lbl));
        
        X_Ind = 1:length(grid_vis.X(:,1));
        Y_Ind = 1:length(grid_vis.Y(1,:));
        %imagesc_sp = imagesc(ax_cur, grid_vis.X(:,1), grid_vis.Y(1,:), dataij');
        %imagesc_sp = imagesc(ax_cur, X_Ind, Y_Ind, dataij');
        %hold on
        ctfsc_sp = contourf(ax_cur, grid_vis.X, grid_vis.Y, dataij);
                
        set(ax_cur, 'xticklabel',{[]});
        set(ax_cur, 'yticklabel',{[]});
        axis(ax_cur, 'xy')
        pbaspect(ax_cur, [1 1 1])
        caxis([-0.1, 0.1])
        colorbar;
        colormap('jet')
        
        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , Fig_outputpath_lbl{lbl}, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
    
    Fig_outputpath_lbl = {'_Uf', '_Vf'};
    for lbl = 1:2
        FigObj = figure('Name', 'Fields_UV', 'NumberTitle','on', 'pos', [200 200 800 800]);
        FigObj_Name = sprintf('U: h= %d days', SamplingInterval_vis)
        ax_cur = gca;
        
        % Happen that 1,2 are also velocity component
        dataij = squeeze(theta_MAP_ij(:,:,lbl));
        
        X_Ind = 1:length(grid_vis.X(:,1));
        Y_Ind = 1:length(grid_vis.Y(1,:));
        %imagesc_sp = imagesc(ax_cur, grid_vis.X(:,1), grid_vis.Y(1,:), dataij');
        %imagesc_sp = imagesc(ax_cur, X_Ind, Y_Ind, dataij');
        %hold on
        
        ctf_Value = linspace(-0.1, 0.1, 20);
        if lbl == 1
            dataij = theta_MAP_ij_f.u;
        elseif lbl == 2
            dataij = theta_MAP_ij_f.v;
        end
        ctfsc_sp = contourf(ax_cur, grid_vis_f.X, grid_vis_f.Y, dataij, ctf_Value);
        
        set(ax_cur, 'xticklabel',{[]});
        set(ax_cur, 'yticklabel',{[]});
        axis(ax_cur, 'xy')
        pbaspect(ax_cur, [1 1 1])
        caxis([-0.1, 0.1])
        colorbar;
        colormap('jet')
        
        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , Fig_outputpath_lbl{lbl}, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
   
    Fig_outputpath_lbl = {'_logK11', '_logK22'};
    for lbl = 1:2
        FigObj = figure('Name', 'Fields_UV', 'NumberTitle','on', 'pos', [200 200 800 800]);
        FigObj_Name = sprintf('U: h= %d days', SamplingInterval_vis)
        ax_cur = gca;
        
        % Happen that 7,8 are also diffusivity component
        dataij = squeeze(log10(theta_MAP_ij(:,:,lbl+6)));
        
        X_Ind = 1:length(grid_vis.X(:,1));
        Y_Ind = 1:length(grid_vis.Y(1,:));
        imagesc_sp = imagesc(ax_cur, X_Ind, Y_Ind, dataij');
        
        %ctfsc_sp = contourf(ax_cur, grid_vis.X, grid_vis.Y, dataij);
        
        set(ax_cur, 'xticklabel',{[]});
        set(ax_cur, 'yticklabel',{[]});
        axis(ax_cur, 'xy')
        pbaspect(ax_cur, [1 1 1])
        caxis([0, 5])
        colorbar;
        colormap('jet')

        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , Fig_outputpath_lbl{lbl}, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
  
   Fig_outputpath_lbl = {'_A11', '_A12', '_A21', '_A22'};
    for lbl = 1:4
        FigObj = figure('Name', 'Fields_UV', 'NumberTitle','on', 'pos', [200 200 800 800]);
        FigObj_Name = sprintf('U: h= %d days', SamplingInterval_vis)
        ax_cur = gca;
        
        % Happen that 7,8 are also diffusivity component
        dataij = squeeze(theta_MAP_ij(:,:,lbl+2));
        
        X_Ind = 1:length(grid_vis.X(:,1));
        Y_Ind = 1:length(grid_vis.Y(1,:));
        imagesc_sp = imagesc(ax_cur, X_Ind, Y_Ind, dataij');
        
        %ctfsc_sp = contourf(ax_cur, grid_vis.X, grid_vis.Y, dataij);
        
        set(ax_cur, 'xticklabel',{[]});
        set(ax_cur, 'yticklabel',{[]});
        axis(ax_cur, 'xy')
        pbaspect(ax_cur, [1 1 1])
        caxis([-1E-7, 1E-7])
        colorbar;
        colormap('jet')

        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , Fig_outputpath_lbl{lbl}, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
    %%
    Fig_outputpath_lbl = {'_detA'};
    for lbl = 1:1
        FigObj = figure('Name', 'Fields_detA', 'NumberTitle','on', 'pos', [200 200 800 800]);
        FigObj_Name = sprintf('U: h= %d days', SamplingInterval_vis)
        ax_cur = gca;
        
        % Happen that 7,8 are also diffusivity component
        dij = theta_MAP_ij ;
        
        dataij = squeeze(dij(:,:,3).*dij(:,:,6)-dij(:,:,4)*dij(:,:,5));
        
        %dataij = sign(dataij);
        
        X_Ind = 1:length(grid_vis.X(:,1));
        Y_Ind = 1:length(grid_vis.Y(1,:));
        imagesc_sp = imagesc(ax_cur, X_Ind, Y_Ind, dataij');
        
        %ctfsc_sp = contourf(ax_cur, grid_vis.X, grid_vis.Y, dataij);
        
        set(ax_cur, 'xticklabel',{[]});
        set(ax_cur, 'yticklabel',{[]});
        axis(ax_cur, 'xy')
        pbaspect(ax_cur, [1 1 1])
        %caxis([0, 5])
        colorbar;
        colormap('jet')

        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , Fig_outputpath_lbl{lbl}, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
end



if (contains(RunProfile, 'QGM2') && contains(RunProfile, 'Fields_vs_SIV'))
    SamplingInterval_vis_List = [1, 16; 64, 128];
    SamplingInterval_vis_List = [160];
    %jet_cmap = colormap('jet');
    
    % ------- Part 1 -------- %
    ax_cur_sp = cell(size(SamplingInterval_vis_List));
    FigObj = figure('Name', 'Fields_Streamline', 'NumberTitle','on', 'pos', [200 200 1200 800]);
    for SIV_ind = 1:length(SamplingInterval_vis_List(:))
        % Create subplot
        ax_cur_sp{SIV_ind} = subplot(size(SamplingInterval_vis_List, 1), size(SamplingInterval_vis_List, 2), SIV_ind);
        ax_cur_cur = gca;
        
        % Load a specfic field
        SamplingInterval_vis = SamplingInterval_vis_List(SIV_ind);
        nTrial = 0;
        
        % Load Inferred Fields
        DataTypeString = 'BayesSampleStat';
        run('Script_Filenames');
        
        BSS_DATA = load(filename_BSS);
        
        % Evaluate the required fields for visualisation
        burnin_bgn = round(BSS_DATA.MCMC_Param.Nsteps_pc *0.50);
        
        theta_samples = BSS_DATA.theta_Stat.theta_store(:, :, burnin_bgn:end);
        %[Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);
        [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_FieldsKPolar(theta_samples);
        % Take squares of  sigma_1, sigma_2 into
        Fields(:, 3:4, :) = Fields(:, 3:4, :).^2;
        
        MeanFieldsMerged = {mean(Fields, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
        
        %  Visualise the Fields and Position of the cell points
        Mesh = BSS_DATA.Mesh_Struct.Mesh;
        RectMesh_Param = BSS_DATA.Mesh_Struct.RectMesh_Param;
        
        grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
        
        % Reorganise the data
        MeanFields_CellAvg = MeanFieldsMerged{1};
        MeanFields_CellAvg_ij = convert_cellk_to_cellij(MeanFields_CellAvg, RectMesh_Param);
        
        U_vis = MeanFields_CellAvg_ij(:,:,1)*1E2;  % convert into cgs
        V_vis = MeanFields_CellAvg_ij(:,:,2)*1E2;  % convert into cgs
        sigma_isq_vis = MeanFields_CellAvg_ij(:,:,3);
        
        UV_mag = sqrt(U_vis.^2+V_vis.^2);
        %U_vis = U_vis./UV_mag;
        %V_vis = V_vis./UV_mag;

        base = 10;
        UV_magorder = floor(log(abs(UV_mag))./log(base));
        
        orderlv_List = [min(UV_magorder(:)):max_cur(UV_magorder(:))];
        orderlv_List = [-5:1];

        sigma_isq_vis = log10(sigma_isq_vis);
        
        %sigma1sq_vis = flipud(sigma1sq_vis);
        %img_str = imagesc(ax_cur_cur, grid_vis.X(:,1), grid_vis.Y(1,:), sigma1sq_vis');
        img_str = imagesc(ax_cur_cur, grid_vis.X(:,1), grid_vis.Y(1,:), UV_mag');
        %img_str = imagesc(grid_vis.X(:,1), grid_vis.Y(1,:), sigma_isq_vis');

        ax_curis(ax_cur_cur, 'xy')
        cax_curis([0, 1.5]);  % Hardcode colourscale
        %colormap('gray')
        %colormap('jet')

        
        hold(ax_cur_cur, 'on');
        
        multilevelMag = 0;
        if multilevelMag == 1
            cyan_vector = [0,1,1]; grey_vector = [0.25,0.25,0.25]; white_vector = [0.85,0.85,0.85];
            colorind_List = linspace(0, 1, length(orderlv_List));
            colorsc_List = zeros(size(colorind_List));
            qv_str_List = cell(size(colorind_List));
            
            for lv = 1:length(orderlv_List)
                colorsc = 1-colorind_List(lv);
                colorsc_List(lv) = colorsc;
                
                maglv = orderlv_List(lv);
                U_vis_norm = U_vis;
                V_vis_norm = V_vis;
                U_vis_norm(UV_magorder ~= maglv) = 0;
                V_vis_norm(UV_magorder ~= maglv) = 0;
                
                qv_str_List{lv} = quiver(ax_cur_cur, grid_vis.X, grid_vis.Y, U_vis_norm, V_vis_norm);
                quiver(grid_vis.X, grid_vis.Y, U_vis_norm, V_vis_norm);
                set(qv_str_List{lv}, 'Color', colorsc*[0.8,0.8,0.8], 'LineStyle','-', 'LineWidth',lv/4)
            end
            
            qv_arr = []; legend_List = {};
            for lv = 1:length(orderlv_List)
                qv_arr = [qv_arr, qv_str_List{lv}];
                legend_List{lv} = sprintf('$M= %i$ ', orderlv_List(lv));
            end
            
            if (length(SamplingInterval_vis_List(:)) == 1)
                title(ax_cur_cur, ['$u_{avg} [cm/s]$ and $\log_{10}(\frac{{\sigma_1}^2} {1 \, m^2/s})$'] );
            end
            %xlabel(ax_cur_cur, sprintf('sampling interval $h=%i$ days', SamplingInterval_vis))
            xlabel(ax_cur_cur, sprintf('$h=%i$ days', SamplingInterval_vis))
            set(ax_cur_cur,'YTickLabel',[], 'XTickLabel',[]);
            
            pbaspect([1 1 1])
            %colorbar;
            run('Script_ax_curesConfig.m')
            
            if (length(SamplingInterval_vis_List(:)) > 1)
                % Use the final subplot to construct the legend and colourbar
                %legend(qv_arr, legend_List, 'Position', [0.15, 0.025, (1-0.15*2), 0.02], 'Orientation', 'vertical')
                legend(qv_arr, legend_List, 'Position', [0.05, 0.5-0.1, 0.05*2, 0.1*2], 'Orientation', 'vertical')
                
                globalcb = colorbar('Position', [0.15, 0.04, (1-0.15*2), 0.02], 'Orientation', 'horizontal', 'FontSize', 14);
                %set(get(globalcb,'title'),'string','$\log_{10}(\frac{{\sigma_1}^2} {1 \, m^2/s})$');
                
                run('Script_ax_curesConfig.m')
            end
        else
            U_vis_norm = U_vis./UV_mag;
            V_vis_norm = V_vis./UV_mag;

            qv_str = quiver_angle(ax_cur_cur, grid_vis.X, grid_vis.Y, U_vis_norm, V_vis_norm);
            
            xlabel(ax_cur_cur, sprintf('$h=%i$ days', SamplingInterval_vis))
            set(ax_cur_cur,'YTickLabel',[], 'XTickLabel',[]);
            
            pbaspect([1 1 1])
            %colorbar;
            run('Script_ax_curesConfig.m')
            
            if SIV_ind == 1
                title('Velocity Field');
                globalcb = colorbar('Position', [0.15, 0.04, (1-0.15*2), 0.02], 'Orientation', 'horizontal', 'FontSize', 12);
                globalcb.Label.Interpreter = 'latex';
            end
        end
    end

    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    if (length(SamplingInterval_vis_List(:)) >= 2)
        Figure_param = ['_diffSIv', '_Idx', num2str(Nx_cell_ARG)];
    else
        Figure_param = ['_SIv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    end
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_UVMag', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    % ------- Part 1b: sigma1,2, phi -------- %
    for sigma_i = 1:2
        ax_cur_sp = cell(size(SamplingInterval_vis_List));
        FigObj = figure('Name', ['Fields_sigma', num2str(sigma_i)], 'NumberTitle','on', 'pos', [200 200 1200 800]);
        for SIV_ind = 1:length(SamplingInterval_vis_List(:))
            % Create subplot
            ax_cur_sp{SIV_ind} = subplot(size(SamplingInterval_vis_List, 1), size(SamplingInterval_vis_List, 2), SIV_ind);
            ax_cur_cur = gca;
            
            % Load a specfic field
            SamplingInterval_vis = SamplingInterval_vis_List(SIV_ind);
            nTrial = 0;
            
            % ------- Part 1 -------- %
            % Load Inferred Fields
            DataTypeString = 'BayesSampleStat';
            run('Script_Filenames');
            
            BSS_DATA = load(filename_BSS);
            
            % Evaluate the required fields for visualisation
            burnin_bgn = round(BSS_DATA.MCMC_Param.Nsteps_pc *0.50);
            
            theta_samples = BSS_DATA.theta_Stat.theta_store(:, :, burnin_bgn:end);
            %[Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);
            [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_FieldsKPolar(theta_samples);
            % Take squares of  sigma_1, sigma_2 into
            Fields(:, 3:4, :) = Fields(:, 3:4, :).^2;
            
            MeanFieldsMerged = {mean(Fields, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
            
            %  Visualise the Fields and Position of the cell points
            Mesh = BSS_DATA.Mesh_Struct.Mesh;
            RectMesh_Param = BSS_DATA.Mesh_Struct.RectMesh_Param;
            
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
            
            % Reorganise the data
            MeanFields_CellAvg = MeanFieldsMerged{1};
            MeanFields_CellAvg_ij = convert_cellk_to_cellij(MeanFields_CellAvg, RectMesh_Param);
            
            U_vis = MeanFields_CellAvg_ij(:,:,1)*1E2;  % convert into cgs
            V_vis = MeanFields_CellAvg_ij(:,:,2)*1E2;  % convert into cgs
            sigma_isq_vis = MeanFields_CellAvg_ij(:,:,3+(sigma_i-1));
            
            sigma_isq_vis = log10(sigma_isq_vis);
            
            %sigma1sq_vis = flipud(sigma1sq_vis);
            img_str = imagesc(ax_cur_cur, grid_vis.X(:,1), grid_vis.Y(1,:), sigma_isq_vis');
            ax_curis(ax_cur_cur, 'xy')
            cax_curis([0, 5]);  % Hardcode colourscale
            %colormap('gray')
            colormap('jet')
            
            hold(ax_cur_cur, 'on');
            
            phi_x = cos(MeanFields_CellAvg_ij(:,:,5) + (sigma_i-1)*pi/2);
            phi_y = sin(MeanFields_CellAvg_ij(:,:,5) + (sigma_i-1)*pi/2);
            
            qv_str = quiver_angle(ax_cur_cur, grid_vis.X, grid_vis.Y, phi_x, phi_y);
            
            xlabel(ax_cur_cur, sprintf('$h=%i$ days', SamplingInterval_vis))
            set(ax_cur_cur,'YTickLabel',[], 'XTickLabel',[]);
            pbaspect(ax_cur_cur, [1 1 1]);
            
            %colorbar;
            if SIV_ind == 1
                sigmai_text = sprintf('{\\sigma_%i}^2', sigma_i);
                title_text = ['$\log_{10} \left( \frac{', sigmai_text, '}{1 \, m^2/s} \right)$'];
                title(ax_cur_cur, title_text);
                globalcb = colorbar(ax_cur_cur, 'Position', [0.15, 0.04, (1-0.15*2), 0.02], 'Orientation', 'horizontal', 'FontSize', 12);
                globalcb.Label.Interpreter = 'latex';
            end
            run('Script_ax_curesConfig.m')
        end
        
        
        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        if (length(SamplingInterval_vis_List(:)) >= 2)
            Figure_param = ['_diffSIv', '_Idx', num2str(Nx_cell_ARG)];
        else
            Figure_param = ['_SIv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
        end
       
        Fig_outputpath = [ Fig_outputfolder , RunProfile , sprintf( '_sigma%isqMag', sigma_i), Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
    
%     %% ------- Part 1c: sigma1,2, phi -------- %
%     ax_cur_sp = cell(size(SamplingInterval_vis_List));
%     FigObj = figure('Name', ['Fields_sigma', num2str(sigma_i)], 'NumberTitle','on', 'pos', [200 200 1200 800]);
%     for SIV_ind = 1:length(SamplingInterval_vis_List(:))
%         % Create subplot
%         ax_cur_sp{SIV_ind} = subplot(size(SamplingInterval_vis_List, 1), size(SamplingInterval_vis_List, 2), SIV_ind);
%         ax_cur_cur = gca;
%         
%         % Load a specfic field
%         SamplingInterval_vis = SamplingInterval_vis_List(SIV_ind);
%         nTrial = 0;
%         
%         % ------- Part 1 -------- %
%         % Load Inferred Fields
%         DataTypeString = 'BayesSampleStat';
%         run('Script_Filenames');
%         
%         BSS_DATA = load(filename_BSS);
%         
%         % Evaluate the required fields for visualisation
%         burnin_bgn = round(BSS_DATA.MCMC_Param.Nsteps_pc *0.50);
%         
%         theta_samples = BSS_DATA.theta_Stat.theta_store(:, :, burnin_bgn:end);
%         %[Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);
%         [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_FieldsKPolar(theta_samples);
%         % Take squares of  sigma_1, sigma_2 into
%         Fields(:, 3:4, :) = Fields(:, 3:4, :).^2;
%         
%         MeanFieldsMerged = {mean(Fields, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
%         
%         %  Visualise the Fields and Position of the cell points
%         Mesh = BSS_DATA.Mesh_Struct.Mesh;
%         RectMesh_Param = BSS_DATA.Mesh_Struct.RectMesh_Param;
%         
%         grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
%         
%         % Reorganise the data
%         MeanFields_CellAvg = MeanFieldsMerged{1};
%         MeanFields_CellAvg_ij = convert_cellk_to_cellij(MeanFields_CellAvg, RectMesh_Param);
%         
%         U_vis = MeanFields_CellAvg_ij(:,:,1)*1E2;  % convert into cgs
%         V_vis = MeanFields_CellAvg_ij(:,:,2)*1E2;  % convert into cgs
%         
%         UV_mag = sqrt(U_vis.^2 + V_vis.^2);
%         sigma_isq_vis = MeanFields_CellAvg_ij(:,:,3);
%         
%         sigma_isq_vis = log10(sigma_isq_vis);
%         
%         %sigma1sq_vis = flipud(sigma1sq_vis);
%         %img_str = imagesc(ax_cur_cur, grid_vis.X(:,1), grid_vis.Y(1,:), sigma1sq_vis');
%         img_str = imagesc(ax_cur_cur, grid_vis.X(:,1), grid_vis.Y(1,:), UV_mag');
%         ax_curis(ax_cur_cur, 'xy')
%         cax_curis([0, 5]);  % Hardcode colourscale
%         %colormap('gray')
%         colormap('jet')
%         
%         hold(ax_cur_cur, 'on');
%         
%         sigma_i = 1;
%         phi_x = cos(MeanFields_CellAvg_ij(:,:,5) + (sigma_i-1)*pi/2);
%         phi_y = sin(MeanFields_CellAvg_ij(:,:,5) + (sigma_i-1)*pi/2);
%         
%         qv_str_phi = quiver_angle(ax_cur_cur, grid_vis.X, grid_vis.Y, phi_x, phi_y);
%         set(qv_str_phi, 'color', 0.5*[1,1,1]')
%         qv_str_vel = quiver_angle(ax_cur_cur, grid_vis.X, grid_vis.Y, U_vis./UV_mag, V_vis./UV_mag);
%         set(qv_str_vel, 'color', 'w')
%         
%         xlabel(ax_cur_cur, sprintf('$h=%i$ days', SamplingInterval_vis))
%         set(ax_cur_cur,'YTickLabel',[], 'XTickLabel',[]);
%         pbaspect(ax_cur_cur, [1 1 1]);
%         
%         %colorbar;
%         if SIV_ind == 1
%             sigmai_text = sprintf('{\\sigma_%i}^2', sigma_i);
%             title_text = ['$\log_{10} \left( \frac{', sigmai_text, '}{1 \, m^2/s} \right)$'];
%             title(ax_cur_cur, title_text);
%             globalcb = colorbar(ax_cur_cur, 'Position', [0.15, 0.04, (1-0.15*2), 0.02], 'Orientation', 'horizontal', 'FontSize', 12);
%             globalcb.Label.Interpreter = 'latex';
%         end
%         run('Script_ax_curesConfig.m')
%     end
%     
%     
%     ImageFormat = 'epsc';
%     Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
%     [status, msg, msgID] = mkdir(Fig_outputfolder);
%     
%     if (length(SamplingInterval_vis_List(:)) >= 2)
%         Figure_param = ['_diffSIv', '_Idx', num2str(Nx_cell_ARG)];
%     else
%         Figure_param = ['_SIv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
%     end
%     
%     Fig_outputpath = [ Fig_outputfolder , RunProfile , sprintf( '_UVphiMag', sigma_i), Figure_param, '.', ImageFormat(1:3)];
%     saveas(FigObj, Fig_outputpath, ImageFormat);
%     
%     
%     
    %%
    % ------- Part 2 -------- %
    ax_cur_sp = cell(size(SamplingInterval_vis_List));
    FigObj = figure('Name', 'RatioRemainNeigh_vis', 'NumberTitle','on', 'pos', [200 200 1200 800]);

    for SIV_ind = 1:length(SamplingInterval_vis_List(:))
        % Create subplot
        ax_cur_sp{SIV_ind} = subplot(size(SamplingInterval_vis_List, 1), size(SamplingInterval_vis_List, 2), SIV_ind);
        ax_cur_cur = gca;
        
        % Load a specfic field
        SamplingInterval_vis = SamplingInterval_vis_List(SIV_ind);
        
        % Load TranDen
        DataTypeString = 'BinnedTrajTransMat';
        run('Script_Filenames');
        BinnedTrajTransMat = load(filename_BTTM);
        
        RatioRemain = BinnedTrajTransMat.RatioRemain;
        RatioRemainNeigh = BinnedTrajTransMat.RatioRemainNeigh;
        
        % Reorganise the data
        RatioRemainNeigh_vis = convert_cellk_to_cellij(RatioRemainNeigh, RectMesh_Param);
        
        
        %sigma1sq_vis = flipud(sigma1sq_vis);
        img_str = imagesc(ax_cur_cur, grid_vis.X(:,1), grid_vis.Y(1,:), RatioRemainNeigh_vis');
        ax_curis(ax_cur_cur, 'xy')
        %colormap('jet')
        
        %if (SIV_ind == length(SamplingInterval_vis_List(:)))
        if (SIV_ind == 2)  %ad-hoc
            hold(ax_cur_cur, 'on');
            [plot_obj, text_obj] = plot_MarkCellCentre(grid_vis, cell_ij_list);
        end
        
        if SIV_ind == 1
            title(ax_cur_cur, ['Remain Proportion'] )
        end
        
        if (length(SamplingInterval_vis_List(:)) == 1)
            xlabel(ax_cur_cur, sprintf('sampling interval $h=%i$ days', SamplingInterval_vis))
        else
            xlabel(ax_cur_cur, sprintf('$h=%i$ days', SamplingInterval_vis))
        end
        set(ax_cur_cur,'YTickLabel',[], 'XTickLabel',[]);
        
        pbaspect([1 1 1])
        colorbar;
        cax_curis([0, 1])
        run('Script_ax_curesConfig.m')
    end
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    if (length(SamplingInterval_vis_List(:)) >= 2)
        Figure_param = ['_diffSIv', '_Idx', num2str(Nx_cell_ARG)];
    else
        Figure_param = ['_SIv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    end
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_RatioRemainNeigh', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end




if (contains(RunProfile, 'QGM2') && contains(RunProfile, 'cmp_vs_SIV'))
    % Main idea:
    % Compare posterior mean of different SIv    
    SamplingInterval_vis_List = [1,2,4,8, 16:16:128];
    
    nTrial_List = [0];
    
    SIv_List_len = length(SamplingInterval_vis_List);
    NTrial = length(nTrial_List);
    
    [Ncells, NVars, MCMC_NSample] = size(theta_store);

    theta_MAP_List = zeros(Ncells, NVars, SIv_List_len);
    RatioRemain_SUM_List = zeros(Ncells, SIv_List_len);
    
    % Load all the data
    for SIv_ind = 1:SIv_List_len
        SamplingInterval_vis = SamplingInterval_vis_List(SIv_ind);
                
        % Compare results across initial conditions
        nTrial = 0;
            
        DataTypeString = 'BayesSampleStat';
        run('Script_Filenames');
        load(filename_BSS);
            
        nTrial_str = num2str(nTrial);
        
        theta_store = theta_Stat.theta_store;
        
        [Ncells, NVars, MCMC_NSample] = size(theta_store);
        
        % Obtain MAP
        logPost_ts = theta_Stat.logPost_ts;
        for cell_k = 1:Ncells
            if mod(log2(cell_k), 4) == 0
                fprintf('Working on %d. \n', cell_k);
            end
            logPost_ts_cell_k = logPost_ts(cell_k, :);
            
            SampInd = find(logPost_ts_cell_k == max(logPost_ts_cell_k(:)), 1);
            
            theta_MAP_cell_k = squeeze(theta_store(cell_k, :, SampInd));
            
            theta_MAP(cell_k, :) = theta_MAP_cell_k;
        end
        
        theta_MAP = convert_theta_U1K0_to_theta_KPolar(theta_MAP);
        [K11, K22, K12] = Kpolar_to_Kcart_vectorised(theta_MAP(:,7,:), theta_MAP(:,8,:), theta_MAP(:,9,:));
        theta_MAP(:,7:9,:) = [K11, K22, K12];
                  
            
        theta_MAP_List(:, :, SIv_ind) = squeeze(theta_MAP);
        
        
        DataTypeString = 'BinnedTrajTransMat';
        run('Script_Filenames');
        load(filename_BTTM);
        
        RatioRemainNeigh(isnan(RatioRemainNeigh)) = 1;
        RatioRemain_SUM_List(:,SIv_ind) = squeeze( RatioRemainNeigh);
        disp(['Finished loading SIv = ', num2str(SamplingInterval_vis)]);
    end
    %%
    SIv_conv_Ncell_ij = zeros(Ncells, 1);   % 5:RatioRemain, RatioRemainNeigh, 
    for cell_k = 1:Ncells
        % 90% of jumps end up in the neigbhouring cells
        SIv_Ind = find((RatioRemain_SUM_List(cell_k, :)>0.9), 1, 'last');
        SIv_conv_Ncell_ij(cell_k) = SamplingInterval_vis_List(SIv_Ind);
    end
    
    %% Visualise data
    NCells_vis = length(cell_k_vis);
    Data_x = SamplingInterval_vis_List;
    
    close all;
    FigObj_List = {};
    for cell_ind = 1:NCells_vis
        cell_k = cell_k_vis(cell_ind);
        cellij = cellk_to_ij(cell_k_vis(cell_ind));
        
        theta_MAP_vs_SIV_cellk = squeeze(theta_MAP_List(cell_k, :, :));
        
                
        FigObj_Name = sprintf('i: %d, j: %d; k: %d', cellij(1), cellij(2), cell_k)
        FigObj_List{cell_ind} = figure('Name', FigObj_Name ,'NumberTitle','on', 'pos', [200 200 800 800]);
        
        FigObj = FigObj_List{cell_ind};
        plotfig = cell(NCells_vis, NVars);

        cell_k = cellij_to_k(cellij(1), cellij(2));
       
        
        Data_Stored = {theta_MAP_vs_SIV_cellk(1, :), theta_MAP_vs_SIV_cellk(2, :), {}; ...
                       theta_MAP_vs_SIV_cellk(4, :), theta_MAP_vs_SIV_cellk(5, :), theta_MAP_vs_SIV_cellk(6, :); ...
                       log10(theta_MAP_vs_SIV_cellk(7, :)), log10(theta_MAP_vs_SIV_cellk(8, :)), theta_MAP_vs_SIV_cellk(9, :); };          
                   
        title_label = {'$u$ [m/s]', '$v$ [m/s]', ''; ...
                       '$\partial_x u$ [s$^{-1}$]', '$\partial_y u$ [s$^{-1}$]', '$\partial_x v$ [s$^{-1}$]',; ...
                       '$log10(\kappa_{11} [m^2/s]) $', '$log10(\kappa_{22} [m^2/s])$', '$\kappa_{12} [m^2/s]$';};
                   
                   
        ylim_label = {[-0.05, 0.05], [-0.05, 0.05], []; ...
                      [-1E-6, 1E-7], [-5E-7, 5E-7], [-5E-7, 5E-7]; ...
                      [0, 12000], [0, 5], [-2500, 2500];}';   
                  
                                     
        Data_Stored = {theta_MAP_vs_SIV_cellk(1, :), theta_MAP_vs_SIV_cellk(2, :); ...
                       log10(theta_MAP_vs_SIV_cellk(7, :)), log10(theta_MAP_vs_SIV_cellk(8, :)); };
                   
        title_label = {'$u$ [m/s]', '$v$ [m/s]'; ...
                       '$log10(\kappa_{11} [m^2/s]) $', '$log10(\kappa_{22} [m^2/s])$';};
                   
        ylim_label = {[-0.025, 0.025], [-0.025, 0.025]; ...
                      [1, 4], [1, 4];};   
                  
        % Transpose to get the right ordering
        Data_Stored = Data_Stored';
        title_label = title_label';
        ylim_label = ylim_label';
        
        ax_suplots = cell(size(Data_Stored, 1), size(Data_Stored, 2));
        for cmp = 1:length(Data_Stored(:))
            if isempty(Data_Stored{cmp}) ~= 1
                ax_suplots{cmp} = subplot(size(Data_Stored, 1), size(Data_Stored, 2), cmp);
                ax_cur = gca;
                hold(ax_cur, 'on');
                
                Data_y = squeeze(Data_Stored{cmp});
                plotfig{cell_ind, cmp} = plot(ax_cur, Data_x, Data_y);
                xlim(ax_cur, [min(Data_x), max(Data_x)]);
                ylim(ax_cur, ylim_label{cmp});
                ylabel(ax_cur, title_label{cmp});
                xlabel(ax_cur, 'Sampling interval $h$ [days]');

                
                %pbaspect(ax_cur, [1.1 1 1])
                set(ax_cur.XLabel, 'FontSize', 22, 'Interpreter', 'latex')
                set(ax_cur.YLabel, 'FontSize', 22, 'Interpreter', 'latex')

                set(ax_cur.XAxis, 'FontSize', 18)
                set(ax_cur.YAxis, 'FontSize', 18)

                %set(ax_cur,'ytick',[])

                %run('Script_ax_curesConfig.m')
            end
        end
       
        
        %suptitle(sprintf('cell (%d, %d)', cellij(1), cellij(2)))
        
        % Plot a verticle line to show locality assumption breaks down
        for cmp = 1:length(Data_Stored(:))
            ax_cur = ax_suplots{cmp};
            
            hold(ax_cur, 'on');
            
            SIv_conv = SIv_conv_Ncell_ij(cell_k);
            plot(ax_cur, [SIv_conv, SIv_conv], ax_cur.YLim, 'k-.', 'LineWidth',2)
        end
        
        
        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_cellk', num2str(cell_k)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , '_cmpVSSIv', Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
        
    end
    
    
end



if contains(RunProfile, 'QGM2_vs_SIv') %|| (adhoc == 1)
    % Main idea:
    % Compare posterior mean of different SIv    
    SamplingInterval_vis_List = [1, 16:16:160];
    %SamplingInterval_vis_List = [160];
    
    nTrial_List = [1:3];
    
    SIv_List_len = length(SamplingInterval_vis_List);
    NTrial = length(nTrial_List);

    MCMC_Param_List = cell(SIv_List_len, 1);
    %theta_Stat_List = cell(SIv_List_len, 1);
    DownSampling_Param_List = cell(SIv_List_len, 1);
    RatioRemain_List = cell(SIv_List_len, 1);
    RatioRemainNeigh_List = cell(SIv_List_len, 1);
    JumpsCFL_List = cell(SIv_List_len, 1);
    
    Fields_mean_List = cell(SIv_List_len, NTrial);
    Fields_vari_List = cell(SIv_List_len, NTrial);
    
    % Load all the data
    for SIv_ind = 1:SIv_List_len
        SamplingInterval_vis = SamplingInterval_vis_List(SIv_ind);
                
        % Compare results across initial conditions
        for nTrial_ind = 1:NTrial
            nTrial = nTrial_List(nTrial_ind);
            
            DataTypeString = 'BayesSampleStat';
            run('Script_Filenames');
            BSS_DATA = load(filename_BSS);
            
            % Load and Process data
            MCMC_Param_List{SIv_ind, 1} = BSS_DATA.MCMC_Param;
            %theta_Stat_List{nTrial} = BSS_DATA.theta_Stat;
            DownSampling_Param_List{SIv_ind, 1} = BSS_DATA.DownSampling_Param;
            
            %legend_label{end+1} = num2str(DownSampling_Param.SamplingInterval);
            
            burnin_bgn = round(BSS_DATA.MCMC_Param.Nsteps_pc *0.50);
            
            theta_samples = BSS_DATA.theta_Stat.theta_store(:, :, burnin_bgn:end);
            %[Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);
            [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_FieldsKPolar(theta_samples);
            % Take squares of  sigma_1, sigma_2 into 
            Fields(:, 3:4, :) = Fields(:, 3:4, :).^2;            
            
            Fields_mean_List{SIv_ind, nTrial_ind} = {mean(Fields, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
            Fields_vari_List{SIv_ind, nTrial_ind} = {var(Fields, 0, 3), var(GradxFields, 0, 3), var(GradyFields, 0, 3)};
        end        
        
        DataTypeString = 'BinnedTrajTransMat';
        run('Script_Filenames');
        BinnedTrajTransMat = load(filename_BTTM);
        
        RatioRemain_List{SIv_ind} = BinnedTrajTransMat.RatioRemain;
        RatioRemainNeigh_List{SIv_ind} = BinnedTrajTransMat.RatioRemainNeigh;
                
        DataTypeString = 'BinnedTrajData';
        run('Script_Filenames');
        BinnedTrajData = load(filename_BTD);
        JumpsCFL_List{SIv_ind} = BinnedTrajData.TrajJumps_DA.JumpsCFL;
        
        
        RectMesh_Param = BinnedTrajData.Mesh_Struct.RectMesh_Param;
        Mesh = BinnedTrajData.Mesh_Struct.Mesh;
                       
        disp(['Finished loading SIv = ', num2str(SamplingInterval_vis)]);
    end
    
    % Plot graphes
    cellij_to_k = RectMesh_Param.cellij_to_k;

    N_cellij = numel(cell_ij_list);

    close all
    subplot_DIM = size(cell_ij_list);

    %% Re-organise data
    
    % Fields against SIv
    Fields_Ncell_ij = zeros(NTrial, SIv_List_len, N_cellij, 5*3);   % 5:u,v, Kxx, Kyy, Kxy, With Grad
    for cellij_iterator = 1:N_cellij
        cell_ij_ind = cell_ij_list{cellij_iterator};
        cellk_ind = cellij_to_k(cell_ij_ind(1), cell_ij_ind(2));
        
        for SIv_ind = 1:SIv_List_len
            for nTrial_ind = 1:NTrial
                FieldsData = Fields_mean_List{SIv_ind, nTrial_ind};
                Fields_Ncell_ij(nTrial_ind, SIv_ind, cellij_iterator, :) = [FieldsData{1}(cellk_ind, :), FieldsData{2}(cellk_ind, :), FieldsData{3}(cellk_ind, :)];
            end            
        end
    end
     
    
    
    % Sort over different trials
    Fields_Ncell_ij = sort(Fields_Ncell_ij, 1);
    
    % Fraction and JumpsCFL against SIv
    FrcJCFL_Ncell_ij = zeros(SIv_List_len, N_cellij, 5);   % 5:RatioRemain, RatioRemainNeigh, 
    for cellij_iterator = 1:N_cellij
        cell_ij_ind = cell_ij_list{cellij_iterator};
        cellk_ind = cellij_to_k(cell_ij_ind(1), cell_ij_ind(2));
        
        for SIv_ind = 1:SIv_List_len
            RatioRemainData = RatioRemain_List{SIv_ind};
            RatioRemainNeighData = RatioRemainNeigh_List{SIv_ind};

            JCFLData = JumpsCFL_List{SIv_ind};
            max_curJCFLData = max_cur(JCFLData')';
            
            FrcJCFL_Ncell_ij(SIv_ind, cellij_iterator, :) = [RatioRemainData(cellk_ind, :), RatioRemainNeighData(cellk_ind, :), JCFLData(cellk_ind, :), max_curJCFLData(cellk_ind, :)];
        end
    end
    
    SIv_conv_Ncell_ij = zeros(N_cellij, 2);   % 5:RatioRemain, RatioRemainNeigh, 
    for cellij_iterator = 1:N_cellij
        % 50% of jumps end up in the same cell
        SIv_conv_Ncell_ij(cellij_iterator, 1) = find((FrcJCFL_Ncell_ij(:, cellij_iterator, 1)>0.5), 1, 'last');
        
        % 80% of jumps end up in the neighbouring cells
        SIv_conv_Ncell_ij(cellij_iterator, 2) = find((FrcJCFL_Ncell_ij(:, cellij_iterator, 2)>0.9), 1, 'last');
    end
    
    %% Visualise data
    Fields_Ncell_ij_vis = Fields_Ncell_ij;
    
    % Unit conversion into cgs
    uv_ind = [1:2, 6:7, 11:12];
    sigma12_ind = uv_ind+2;
    Fields_Ncell_ij_vis(:,:,:, uv_ind) = Fields_Ncell_ij_vis(:,:,:, uv_ind)*1E2;  % into cgs
    Fields_Ncell_ij_vis(:,:,:, sigma12_ind) = Fields_Ncell_ij_vis(:,:,:, sigma12_ind)*1;  % remains SI
    
    close all
    
    NField = 5;
    if (NField == 9)
        FieldLabel = {'$b_1$', '$b_2$', '${\sigma_{1}}^2$', '${\sigma_{2}}^2$', '$\phi$', 'u_x',  'u_y',  'v_x',  'v_y'};
        Field_ind_List = [1:5, 6:7, 11:12];
    elseif (NField == 5)
        FieldLabel = {'$b_1$ [$cm/s$]', '$b_2$  [$cm/s$]', '${\sigma_{1}}^2$ [$m^2/s$]', '${\sigma_{2}}^2$ [$m^2/s$]', '$\phi [$rad$]$'};
        FieldLabelText = {'u', 'v', 'sigma1', 'sigma2', 'phi'};
        GlobalYlim = {'u', 'v', 'sigma1', 'sigma2', 'phi'};
        
        Field_ind_List = [1:5];
    end
    
    for Field_ind_iterator = 1:length(Field_ind_List)
        Field_ind = Field_ind_List(Field_ind_iterator);
        FieldLabel_cur = FieldLabel{Field_ind_iterator};
        FieldLabelText_cur = FieldLabelText{Field_ind_iterator};
        
        FigObj = figure('Name', FieldLabelText_cur, 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);  
        
        Data_ij_vis_List = cell(size(cell_ij_list));
        YLabelString_List = cell(1, size(cell_ij_list, 2));
        TitleString_List = cell(size(cell_ij_list, 1), 2);

        ind = 0;
        for i =  1:size(Data_ij_vis_List, 1)
            for j = 1:size(Data_ij_vis_List, 2)
                ind = ind + 1;
                Data_ij_vis_List{ind} = squeeze(Fields_Ncell_ij_vis(:, :, ind, Field_ind));
            end
        end
        
        for i =  1:size(cell_ij_list, 1)
            YLabelString_List{i} = ['$j:$ ', num2str( cell_ij_list{i, 1}(2) )];
        end
        for j =  1:size(cell_ij_list, 2)
            TitleString_List{j} = ['$j:$ ', num2str( cell_ij_list{1, j}(1) )];
        end
       
        ax_cur_sp = plot_RClabel(SamplingInterval_vis_List, Data_ij_vis_List, YLabelString_List, TitleString_List);
        
        for cellij_iterator = 1:N_cellij
            ax_cur_spc = ax_cur_sp{cellij_iterator};

            hold(ax_cur_spc, 'on');
%             SIv_conv_ind = SIv_conv_Ncell_ij(cellij_iterator, 1);
%             SIv_conv = SamplingInterval_vis_List(SIv_conv_ind);
%             plot(ax_cur_spc, [SIv_conv, SIv_conv], ax_cur_spc.YLim, 'b-.') 
            
            SIv_conv_ind = SIv_conv_Ncell_ij(cellij_iterator, 2);
            SIv_conv = SamplingInterval_vis_List(SIv_conv_ind);
            plot(ax_cur_spc, [SIv_conv, SIv_conv], ax_cur_spc.YLim, 'k-.', 'LineWidth',2)
        end
        
        % Append title to the first row for the field name
        ax_cur_spc = ax_cur_sp{1, round(size(cell_ij_list, 2)/2)};
        title(ax_cur_spc, {FieldLabel_cur ;ax_cur_spc.Title.String});
        % Append title to the end row for the x label
        ax_cur_spc = ax_cur_sp{end, round(size(cell_ij_list, 2)/2)};
        xlabel(ax_cur_spc, 'sampling interval h (in days)');
        
        ImageFormat = 'epsc';
        Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_vsSIv', '_Idx', num2str(Nx_cell_ARG)]
        Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_cells_', FieldLabelText_cur, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
    FigObj = figure('Name', 'Fraction and JumpsCFL', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    for cellij_iterator = 1:N_cellij
        subplot(subplot_DIM(1), subplot_DIM(2), cellij_iterator);
        Data_plot = squeeze(FrcJCFL_Ncell_ij(:, cellij_iterator, [1,2, 5]));
        
        %boxplot(Data_plot) %, 'Notch','on', 'Labels', string(SamplingInterval_vis_List));
        plot(SamplingInterval_vis_List, Data_plot)
        xlim(SamplingInterval_vis_List([1, end]))
%        legend({'RR','RRN','JCFLmax_cur'});
%        legend({'RR','RRN','JCFLx','JCFLy','JCFLmax_cur'});

        hold on
        % Plot a line of 1
        plot(SamplingInterval_vis_List, ones(size(SamplingInterval_vis_List)), 'k:');

    end
    
    
end


if contains(RunProfile, 'QGM2_vs_nTrial') %|| (adhoc == 1)
    % Main idea:
    % Compare posterior mean of different trial of the same configuration    
    SamplingInterval_vis_List = [1, 16:16:160];

    SIv_List_len = length(SamplingInterval_vis_List);
    
    Max_curRelMeanDiff_ij_List = cell(SIv_List_len, 1);
    MeanFieldsVarMerged_ij_List = cell(SIv_List_len, 1);
    MeanFieldsMerged_ij_List = cell(SIv_List_len, 1);
    Max_curRelMeanDiffKPolar_ij_List = cell(SIv_List_len, 1);
    VarFieldsKPolarMerged_ij_List = cell(SIv_List_len, 1);
    MeanFieldsKPolarMerged_ij_List = cell(SIv_List_len, 1);
    Max_curAbsDiff_ij_List = cell(SIv_List_len, 1);
    Max_curAbsDiffKPolar_ij_List = cell(SIv_List_len, 1);
    
    % Choose SIv
    for SIv_iterator = 1:SIv_List_len
    SamplingInterval_vis = SamplingInterval_vis_List(SIv_iterator);
    
    %% Compare results across initial conditions
    Fields_mean_nTrialList = {};
    Fields_vari_nTrialList = {};
    FieldsKPolar_mean_nTrialList = {};
    FieldsKPolar_vari_nTrialList = {};

    % Compare results across initial conditions    
    MCMC_Param_List = {};
    theta_Stat_List = {};
    DownSampling_Param_List = {};
    
    BSS_DATA_List = cell(length(nTrial_List), 1);
    
    for nTrial = nTrial_List
        DataTypeString = 'BayesSampleStat';
        run('Script_Filenames');
        BSS_DATA = load(filename_BSS);
        BSS_DATA_List{nTrial} = load(filename_BSS);
        
        tic;
        % Load and Process data
        MCMC_Param_List{nTrial} = BSS_DATA.MCMC_Param;
        theta_Stat_List{nTrial} = BSS_DATA.theta_Stat;
        DownSampling_Param_List{nTrial} = BSS_DATA.DownSampling_Param;
        
        %legend_label{end+1} = num2str(DownSampling_Param.SamplingInterval);
        
        burnin_bgn = round(BSS_DATA.MCMC_Param.Nsteps_pc *0.50);
        
        theta_samples = BSS_DATA.theta_Stat.theta_store(:, :, burnin_bgn:end);
        [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);
        [FieldsKPolar, GradxFieldsKPolar, GradyFieldsKPolar] = convert_theta_U1K0_to_FieldsKPolar(theta_samples);

        Fields_mean_nTrialList{end+1} = {mean(Fields, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
        Fields_vari_nTrialList{end+1} = {var(Fields, 0, 3), var(GradxFields, 0, 3), var(GradyFields, 0, 3)};
        
        FieldsKPolar_mean_nTrialList{end+1} = {mean(FieldsKPolar, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
        FieldsKPolar_vari_nTrialList{end+1} = {var(FieldsKPolar, 0, 3), var(GradxFields, 0, 3), var(GradyFields, 0, 3)};
        toc;
    end
     
    %% Compute max_curimal difference
    
    % Merge the fields
    NTrial = length(nTrial_List);
    FieldsMerged = cell(NTrial, 1);
    FieldsKPolarMerged = cell(NTrial, 1);

    FieldsVarMerged = cell(NTrial, 1);
    FieldsKPolarVarMerged = cell(NTrial, 1);
    
    for expt = 1:NTrial
        FieldsMerged{expt} = [Fields_mean_nTrialList{expt}{1}, Fields_mean_nTrialList{expt}{2}, Fields_mean_nTrialList{expt}{3}];
        FieldsVarMerged{expt} = [Fields_vari_nTrialList{expt}{1}, Fields_vari_nTrialList{expt}{2}, Fields_vari_nTrialList{expt}{3}];
        
        FieldsKPolarMerged{expt} = [FieldsKPolar_mean_nTrialList{expt}{1}, FieldsKPolar_mean_nTrialList{expt}{2}, FieldsKPolar_mean_nTrialList{expt}{3}];
        FieldsKPolarVarMerged{expt} = [FieldsKPolar_vari_nTrialList{expt}{1}, FieldsKPolar_vari_nTrialList{expt}{2}, FieldsKPolar_vari_nTrialList{expt}{3}];
    end
    
    % Mean of three mean fields
    MeanFieldsMerged = zeros(size(FieldsMerged{1}));
    MeanFieldsVarMerged = zeros(size(FieldsMerged{1}));
    MeanFieldsKPolarMerged  = zeros(size(FieldsMerged{1}));
    VarFieldsKPolarMerged  = zeros(size(FieldsMerged{1}));

    for expt = 1:NTrial
        MeanFieldsMerged = MeanFieldsMerged + FieldsMerged{expt};
        MeanFieldsVarMerged = MeanFieldsVarMerged + FieldsVarMerged{expt};

        MeanFieldsKPolarMerged = MeanFieldsKPolarMerged + FieldsKPolarMerged{expt};
        VarFieldsKPolarMerged = VarFieldsKPolarMerged + FieldsKPolarVarMerged{expt};
    end
    MeanFieldsMerged = MeanFieldsMerged/NTrial;
    MeanFieldsVarMerged = MeanFieldsVarMerged/NTrial;
    MeanFieldsKPolarMerged = MeanFieldsKPolarMerged/NTrial;
    VarFieldsKPolarMerged = VarFieldsKPolarMerged/NTrial;

    
    % Absolute difference among inferred fields
    Max_curAbsDiff = zeros(size(FieldsMerged{1}));
    Max_curAbsDiffKPolar = zeros(size(FieldsMerged{1}));
    for i = 1:(NTrial-1)
        for j = (i+1):NTrial
            Max_curAbsDiff = max_cur(Max_curAbsDiff, abs(FieldsMerged{j}-FieldsMerged{i}));
            Max_curAbsDiffKPolar = max_cur(Max_curAbsDiffKPolar, abs(FieldsKPolarMerged{j}-FieldsKPolarMerged{i}));
        end
    end
    
    % Difference from the mean inferred fields
    Max_curRelDiff = zeros(size(FieldsMerged{1}));
    for i = 1:NTrial
        Max_curRelDiff = max_cur(Max_curRelDiff, abs(FieldsMerged{i} - MeanFieldsMerged));
    end
    
    % Difference from the mean inferred fields
    Max_curRelDiffKPolar = zeros(size(FieldsKPolarMerged{1}));
    for i = 1:NTrial
        Max_curRelDiffKPolar = max_cur(Max_curRelDiffKPolar, abs(FieldsKPolarMerged{i} - MeanFieldsKPolarMerged));
    end
    
    
    % Assume same mesh for all trials
    Mesh = BSS_DATA.Mesh_Struct.Mesh;
    RectMesh_Param = BSS_DATA.Mesh_Struct.RectMesh_Param;

    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
    
    MeanFieldsMerged_ij = convert_cellk_to_cellij(MeanFieldsMerged, RectMesh_Param);
    MeanFieldsVarMerged_ij = convert_cellk_to_cellij(MeanFieldsVarMerged, RectMesh_Param);
    Max_curAbsDiff_ij = convert_cellk_to_cellij(Max_curAbsDiff, RectMesh_Param);
    Max_curRelDiff_ij = convert_cellk_to_cellij(Max_curRelDiff, RectMesh_Param);
    Max_curRelMeanDiff_ij = Max_curRelDiff_ij./abs(MeanFieldsMerged_ij);
    Max_curAbsMeanDiff_ij = Max_curAbsDiff_ij./abs(MeanFieldsMerged_ij);
    
    % KPolar
    Max_curAbsDiffKPolar_ij = convert_cellk_to_cellij(Max_curAbsDiffKPolar, RectMesh_Param);

    MeanFieldsKPolarMerged_ij = convert_cellk_to_cellij(MeanFieldsKPolarMerged, RectMesh_Param);
    VarFieldsKPolarMerged_ij = convert_cellk_to_cellij(VarFieldsKPolarMerged, RectMesh_Param);
    Max_curRelDiffKPolar_ij = convert_cellk_to_cellij(Max_curRelDiffKPolar, RectMesh_Param);
    Max_curRelMeanDiffKPolar_ij = Max_curRelDiffKPolar_ij./abs(MeanFieldsKPolarMerged_ij);
    
    close all
    plotGraphs = 0;
    if plotGraphs == 1
        for expt = 1:NTrial
            FieldsMerged_ij = convert_cellk_to_cellij(FieldsMerged{expt}, RectMesh_Param);
            FigObj = figure('Name', ['All trials: ', num2str(expt), 'FieldsMerged'], 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
            ax_cur_sp = plot_dataij_k( grid_vis, FieldsMerged_ij(:, :, [1:7, 11:12]) );
        end
        
        for expt = 1:NTrial
            FieldsVarMerged_ij = convert_cellk_to_cellij(FieldsVarMerged{expt}, RectMesh_Param);
            FigObj = figure('Name', ['All trials: ', num2str(expt), 'MeanFieldsVarMerged'], 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
            ax_cur_sp = plot_dataij_k( grid_vis, max_cur(log10(FieldsVarMerged_ij(:, :, [1:7, 11:12])), 0) );
        end
        
        for expt = 1:NTrial
            FieldsKPolarMerged_ij = convert_cellk_to_cellij(FieldsKPolarMerged{expt}, RectMesh_Param);
            FigObj = figure('Name', ['All trials: ', num2str(expt), 'FieldsKPolarMerged'], 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
            ax_cur_sp = plot_dataij_k( grid_vis, FieldsKPolarMerged_ij(:, :, [1:7, 11:12]) );
        end
        
        for expt = 1:NTrial
            FieldsKPolarVarMerged_ij = convert_cellk_to_cellij(FieldsKPolarVarMerged{expt}, RectMesh_Param);
            FigObj = figure('Name', ['All trials: ', num2str(expt), 'MeanFieldsVarMerged'], 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
            ax_cur_sp = plot_dataij_k( grid_vis, max_cur(log10(FieldsKPolarVarMerged_ij(:, :, [1:7, 11:12])), 0) );
        end
        
        FigObj = figure('Name', 'All trials: MeanFieldsMerged', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_cur_sp = plot_dataij_k( grid_vis, MeanFieldsMerged_ij(:, :, [1:7, 11:12]) );
        
        FigObj = figure('Name', 'All trials: MeanFieldsMerged', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_cur_sp = plot_dataij_k( grid_vis, MeanFieldsVarMerged_ij(:, :, [1:7, 11:12]) );
        
        FigObj = figure('Name', 'All trials: Max_curAbsDiff', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_cur_sp = plot_dataij_k( grid_vis, Max_curAbsDiff_ij(:, :, [1:7, 11:12]) );
        
        FigObj = figure('Name', 'All trials: Max_curRelDiff', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_cur_sp = plot_dataij_k( grid_vis, Max_curRelDiff_ij(:, :, [1:7, 11:12]) );
        
        % Output Max_curimal Relative Difference, normalised by the mean
        FigObj = figure('Name', 'All trials: Max_curRelMeanDiff_ij', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_cur_sp = plot_dataij_k( grid_vis, Max_curRelMeanDiff_ij(:, :, [1:7, 11:12]) );
        for k = 1:length(ax_cur_sp)
            ax_cur_sp(k).CLim = [0, 0.3];
        end        
        
        % Output Max_curimal Absolute Difference, normalised by the mean
        FigObj = figure('Name', 'All trials: Max_curAbsMeanDiff_ij', 'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_cur_sp = plot_dataij_k( grid_vis, Max_curAbsMeanDiff_ij(:, :, [1:7, 11:12]) );
        for k = 1:length(ax_cur_sp)
            ax_cur_sp(k).CLim = [0, 0.3];
        end
        % Save the figure
        ImageFormat = 'png';
        Fig_outputfolder = [DataFolder_Figure, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_SIv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)]
        Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_Max_curAbsMeanDiff', Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
    
    % Save for further output
    MeanFieldsMerged_ij_List{SIv_iterator} = MeanFieldsMerged_ij;
    Max_curRelMeanDiff_ij_List{SIv_iterator} = Max_curRelMeanDiff_ij;
    MeanFieldsVarMerged_ij_List{SIv_iterator} = MeanFieldsVarMerged_ij;
    
    Max_curAbsDiff_ij_List{SIv_iterator} = Max_curAbsDiff_ij;
    Max_curAbsDiffKPolar_ij_List{SIv_iterator} = Max_curAbsDiffKPolar_ij;

    MeanFieldsKPolarMerged_ij_List{SIv_iterator} = MeanFieldsKPolarMerged_ij;
    Max_curRelMeanDiffKPolar_ij_List{SIv_iterator} = Max_curRelMeanDiffKPolar_ij;
    VarFieldsKPolarMerged_ij_List{SIv_iterator} = VarFieldsKPolarMerged_ij;
    
    end
   
    %% For publication
    % Target SIv to visualise
    SamplingInterval_vis_List_Target = [1, 16, 64, 128];
    
    SIv_Ind_List_Target = zeros(size(SamplingInterval_vis_List_Target));
    for ind = 1:length(SamplingInterval_vis_List_Target)
        SIv_target = SamplingInterval_vis_List_Target(ind);
        SIv_Ind = find(SamplingInterval_vis_List == SIv_target);
        
        SIv_Ind_List_Target(ind) = SIv_Ind;
    end
    
    TitleString_List = cell(1, length(SIv_Ind_List_Target));
    for c_ind = 1:length(SIv_Ind_List_Target)
        SIv_Ind_target = SIv_Ind_List_Target(c_ind);
        SIv_target = SamplingInterval_vis_List(SIv_Ind_target);
        
        TitleString_List{c_ind} = ['$h = ', num2str(SIv_target), '$ days'];
    end
    
    close all
    %% Visualise MeanFieldsKPolarMerged_ij_List
    % ------ ---------------------------------------------- ------ %
    % ------ Part 1: MeanFieldsKPolarMerged_ij_List ------ %
    
    FigObj = figure('Name', 'All trials: MeanFieldsKPolarMerged_ij_List', 'NumberTitle','on', 'pos', [200 200 1200 800]);
    NFields = 5;    SIv_len = length(SIv_Ind_List_Target);
    Data_ij_vis_List = cell(NFields, SIv_len);

    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
    
    if (NFields == 5)
        NFields_List = 1:5;
        
        row_CLim = {[-4, 1], [-4, 1], [-4, 2], [-4, 2], [-4, 2]}; %, [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 7)
        NFields_List = 1:7;
        
        row_CLim = {[-2, 1], [-2, 1], [-4, 1], [-4, 1], [-2, 1], [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 2)
        NFields_List = 3:4;
        
        row_CLim = {[-4, 1], [-4, 1]};
        YLabelString_List = {'${\sigma_1}^2$', '${\sigma_2}^2$'};
    elseif (NFields == 4)
        NFields_List = 1:4;
        
        row_CLim = {[-4, -1], [-4, -1], [-3, 0], [-3, 0]};
        YLabelString_List = {'$b_1 \, [cm/s]$', '$b_2 \, [cm/s]$', '${\sigma_1}^2 \, [m^2/s]$', '${\sigma_2}^2 \, [m^2/s]$'};
    end
    
    % Read in Data
    for Fld_it = 1:length(NFields_List)
        Fld_ind = NFields_List(Fld_it);
        for SIv_Ind = 1:length(SIv_Ind_List_Target)
            if Fld_ind <= 2
                % Convert into cgs
                Data_ij_vis_List{Fld_it, SIv_Ind} = (100*(squeeze(MeanFieldsKPolarMerged_ij_List{SIv_Ind}(:, :, Fld_ind)))); % sigma_1, sigma_2, phi
            else
                Data_ij_vis_List{Fld_it, SIv_Ind} = (squeeze(MeanFieldsKPolarMerged_ij_List{SIv_Ind}(:, :, Fld_ind))); % sigma_1, sigma_2, phi
            end
        end
    end

    %ax_cur_sp = pcolor_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    ax_cur_sp = imagesc_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    
    if NFields >= 5
        % Admend the colormap for cyclic scale
        for c_ind = 1:length(SIv_Ind_List_Target)
            colormap(ax_cur_sp{5, c_ind}, hsv);
        end
    end

    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_vsSIv', '_Idx', num2str(Nx_cell_ARG)]
    Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_MeanFieldsUVKPolar', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    % ------ ---------------------------------------------- ------ %
    % ------ Part 2: Visualise Max_curRelMeanDiffKPolar_ij_List ------ %
    
    FigObj = figure('Name', 'All trials: Max_curRelMeanDiffKPolar_ij_List', 'NumberTitle','on', 'pos', [200 200 1200 800]);
    NFields = 4;    
    SIv_len = length(SIv_Ind_List_Target);
    Data_ij_vis_List = cell(NFields, SIv_len);
    
    if (NFields == 5)
        NFields_List = 1:5;
        
        row_CLim = {[-2, 1], [-2, 1], [-4, 1], [-4, 1], [-2, 1]}; %, [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'}; %, 'dudx', 'dudy'};
    elseif (NFields == 7)
        NFields_List = 1:7;
        
        row_CLim = {[-2, 1], [-2, 1], [-4, 1], [-4, 1], [-2, 1], [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 2)
        NFields_List = 3:4;
        
        row_CLim = {[-4, 1], [-4, 1]};
        YLabelString_List = {'${\sigma_1}^2$', '${\sigma_2}^2$'};
    elseif (NFields == 4)
        NFields_List = 1:4;
        
        row_CLim = {[-2, 1], [-2, 1], [-4, 1], [-4, 1]};
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$'};    
    end
    
    % Read in Data
    for Fld_it = 1:length(NFields_List)
        Fld_ind = NFields_List(Fld_it);
        for SIv_Ind = 1:length(SIv_Ind_List_Target)
            Data_ij_vis_List{Fld_it, SIv_Ind} = log10(squeeze(Max_curRelMeanDiffKPolar_ij_List{SIv_Ind}(:, :, Fld_ind))); % sigma_1, sigma_2, phi
        end
    end

    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
    

    %ax_cur_sp = pcolor_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    ax_cur_sp = imagesc_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    
    % Admend the colormap for cyclic scale for phi
    if (NFields >= 5)
        for c_ind = 1:length(SIv_Ind_List_Target)
            colormap(ax_cur_sp{5, c_ind}, hsv);
        end
    end
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_vsSIv', '_Idx', num2str(Nx_cell_ARG)]
    Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_Max_curRelMeanDiffKPolar', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
    % ------ ---------------------------------------------- ------ %
    %% ------ Part 3: Visualise VarFieldsKPolarMerged_ij_List ------ %
    
    FigObj = figure('Name', 'All trials: VarFieldsKPolarMerged', 'NumberTitle','on', 'pos', [200 200 1200 800]);
    NFields = 4;    SIv_len = length(SIv_Ind_List_Target);
    Data_ij_vis_List = cell(NFields, SIv_len);

    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
    
    if (NFields == 5)
        NFields_List = 1:5;
        
        row_CLim = {[-4, 1], [-4, 1], [-4, 2], [-4, 2], [-4, 2]}; %, [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 7)
        NFields_List = 1:7;
        
        row_CLim = {[-2, 1], [-2, 1], [-4, 1], [-4, 1], [-2, 1], [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 2)
        NFields_List = 3:4;
        
        row_CLim = {[-4, 1], [-4, 1]};
        YLabelString_List = {'${\sigma_1}^2$', '${\sigma_2}^2$'};
    elseif (NFields == 4)
        NFields_List = 1:4;
        
        row_CLim = {[-4, -1], [-4, -1], [-3, 0], [-3, 0]};
        YLabelString_List = {'$b_1 \, [cm/s]$', '$b_2 \, [cm/s]$', '${\sigma_1}^2 \, [m^2/s]$', '${\sigma_2}^2 \, [m^2/s]$'};
    end
    
    % Read in Data
    for Fld_it = 1:length(NFields_List)
        Fld_ind = NFields_List(Fld_it);
        for SIv_Ind = 1:length(SIv_Ind_List_Target)
            if Fld_ind <= 2
                % Convert into cgs
                Data_ij_vis_List{Fld_it, SIv_Ind} = log10(100*sqrt(squeeze(VarFieldsKPolarMerged_ij_List{SIv_Ind}(:, :, Fld_ind)))); % sigma_1, sigma_2, phi
            else
                Data_ij_vis_List{Fld_it, SIv_Ind} = log10(sqrt(squeeze(VarFieldsKPolarMerged_ij_List{SIv_Ind}(:, :, Fld_ind)))); % sigma_1, sigma_2, phi
                
            end
        end
    end

    %ax_cur_sp = pcolor_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    ax_cur_sp = imagesc_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    
    if NFields >= 5
        % Admend the colormap for cyclic scale
        for c_ind = 1:length(SIv_Ind_List_Target)
            colormap(ax_cur_sp{5, c_ind}, hsv);
        end
    end
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_vsSIv', '_Idx', num2str(Nx_cell_ARG)]
    Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_VarFieldsKPolar', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
   
    
    % ------ ---------------------------------------------- ------ %
    %% ------ Part 4: Visualise Max_curDiffKPolar_ij_List ------ %
    FigObj = figure('Name', 'All trials: Max_curAbsDiffKPolar_ij_List', 'NumberTitle','on', 'pos', [200 200 1200 800]);
    NFields = 5;    
    SIv_len = length(SIv_Ind_List_Target);
    Data_ij_vis_List = cell(NFields, SIv_len);
    
    if (NFields == 5)
        NFields_List = 1:5;
        
        row_CLim = {[-4, 1], [-4, 1], [-4, 2], [-4, 2], [-4, 2]}; %, [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 7)
        NFields_List = 1:7;
        
        row_CLim = {[-2, 1], [-2, 1], [-4, 1], [-4, 1], [-2, 1], [-2, 1], [-2, 1] };
        YLabelString_List = {'$b_1$', '$b_2$', '${\sigma_1}^2$', '${\sigma_2}^2$', '$\phi$'};
    elseif (NFields == 2)
        NFields_List = 3:4;
        
        row_CLim = {[-4, 1], [-4, 1]};
        YLabelString_List = {'${\sigma_1}^2$', '${\sigma_2}^2$'};
    end
    
    % Read in Data
    for Fld_it = 1:length(NFields_List)
        Fld_ind = NFields_List(Fld_it);
        for SIv_Ind = 1:length(SIv_Ind_List_Target)
            Data_ij_vis_List{Fld_it, SIv_Ind} = log10(squeeze(Max_curAbsDiffKPolar_ij_List{SIv_Ind}(:, :, Fld_ind))); % sigma_1, sigma_2, phi
            Data_ij_vis_List{Fld_it, SIv_Ind} = (squeeze(Max_curAbsDiffKPolar_ij_List{SIv_Ind}(:, :, Fld_ind))); % sigma_1, sigma_2, phi
        end
    end

    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
    

    %ax_cur_sp = pcolor_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    ax_cur_sp = imagesc_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);
    
    % Admend the colormap for cyclic scale for phi
    if (NFields >= 5)
        for c_ind = 1:length(SIv_Ind_List_Target)
            colormap(ax_cur_sp{5, c_ind}, hsv);
        end
    end
    
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_vsSIv', '_Idx', num2str(Nx_cell_ARG)]
    Fig_outputpath = [ Fig_outputfolder , RunProfile ,'_Max_curAbsDiff', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
        
    % Output to screen
    for SIv_Ind = 1:length(SIv_Ind_List_Target)
        SIv = SamplingInterval_vis_List_Target(SIv_Ind);
        disp(['Sampling Interval: ', num2str(SIv)]);
        
        for cmp = 1:size(Data_ij_vis_List)
            SIv = SamplingInterval_vis_List_Target(SIv_Ind);
            disptxt = ['Component ', num2str(cmp), ': ',  num2str(mean(Data_ij_vis_List{cmp, SIv_Ind}(:)))];
            disp(disptxt);
        end
    end
end


if strcmp(RunProfile, 'cosine_TC_pwc') || strcmp(RunProfile, 'cosine_TC')
    %% Output plot of inferred quantities vs Sampling Intervals
    
    Nexpt = length(theta_Stat_List);
    
    Fields_mean_nTrialList = {};
    Fields_vari_nTrialList = {};
        
    for expt = 1:Nexpt
        tic;
        % Load and Process data
        MCMC_Param = MCMC_Param_List{expt};
        theta_Stat = theta_Stat_List{expt};
        DownSampling_Param = DownSampling_Param_List{expt};
        
        %legend_label{end+1} = num2str(DownSampling_Param.SamplingInterval);
        
        burnin_bgn = round(MCMC_Param.Nsteps_pc *0.25);
        
        theta_samples = theta_Stat.theta_store(:, :, burnin_bgn:end);
        [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);

        Fields_mean_nTrialList{end+1} = {mean(Fields, 3), mean(GradxFields, 3), mean(GradyFields, 3)};
        Fields_vari_nTrialList{end+1} = {var(Fields, 0, 3), var(GradxFields, 0, 3), var(GradyFields, 0, 3)};
        toc;
    end
          
    % Plot against y only
    % Assumed same grid size
    grid_y_bdy = Mesh_Struct.RectMesh_Param.grid_y';
    y_coord = 0.5*(grid_y_bdy(1:end-1) + grid_y_bdy(2:end));
    
    % Compute exact quantities
    U_exact = veloc_fldStruct.u(0*y_coord', y_coord');
    kappa_exact = zeros(3, length(y_coord));  % 3: for KXX, KYY, KXY
    for j_y = 1:length(y_coord)
        kappa_yj = kappa_fldStruct.kappa(0, y_coord(j_y));
        kappa_exact(1, j_y)  = kappa_yj(1,1);
        kappa_exact(2, j_y)  = kappa_yj(2,2);
        kappa_exact(3, j_y)  = kappa_yj(1,2);
    end
    
    % Compute left-edge and right-edge field
    Fields_yl_List = cell(size(Fields_mean_nTrialList));
    Fields_yr_List = cell(size(Fields_mean_nTrialList));
    
    grid_y_bdy_l = grid_y_bdy(1:end-1);
    grid_y_bdy_r = grid_y_bdy(2:end);
    
    diffy_halfl  = (y_coord - grid_y_bdy_l);
    diffy_halfr  = (grid_y_bdy_r - y_coord);
    
    for k = 1:length(Fields_mean_nTrialList)
        Fields_mean = Fields_mean_nTrialList{k}{1};
        GradxFields_mean = Fields_mean_nTrialList{k}{2};
        GradyFields_mean = Fields_mean_nTrialList{k}{3};
        
        Fields_yl_List{k} = Fields_mean - diffy_halfl.*GradyFields_mean;
        Fields_yr_List{k} = Fields_mean + diffy_halfr.*GradyFields_mean;
    end
       
    
    % Plot fields under different sampling interval
    %expt_plot = [1, 18, 34, 50, 66, 66+32];
    expt_plot = [1, 34:32:(34+3*32)];
    expt_plot = [1, 66:64:(34+3*32)];

    Nexpt_plot = length(expt_plot);
        
    % Output screen
    FigObj = figure('Name', 'vs SIv', 'NumberTitle', 'off', 'pos', [200 200 1200 800]);
    FigObjax_cur = ax_cures;
    
    plot_ind = 0;
    for expt = expt_plot
        SamplingInterval = DownSampling_Param_List{expt}.SamplingInterval;
        
        plot_ind = plot_ind + 1;
        subplotax_cur = subplot(Nexpt_plot,2, plot_ind);
        plot(subplotax_cur, ...
            [grid_y_bdy_l, grid_y_bdy_r]', [Fields_yl_List{expt}(:, 1), Fields_yr_List{expt}(:, 1)]', 'b', ...
            [grid_y_bdy_l, grid_y_bdy_r]', [Fields_yl_List{expt}(:, 2), Fields_yr_List{expt}(:, 2)]', 'r')
        hold on
        plot(subplotax_cur, y_coord, U_exact(1,:), 'b:', ...
            y_coord, U_exact(2,:), 'r:')
        xlabel('$y$');
        %legend('u', 'v', 'u_e', 'v_e')
        ylabel(['h = ', num2str(SamplingInterval)]);
        ylim([-1.2, 1.2])
        run('Script_ax_curesConfig.m')

        
        plot_ind = plot_ind + 1;
        subplotax_cur = subplot(Nexpt_plot,2, plot_ind);
        plot(subplotax_cur, ...
            [grid_y_bdy_l, grid_y_bdy_r]', [Fields_yl_List{expt}(:, 3), Fields_yr_List{expt}(:, 3)]', 'b', ...
            [grid_y_bdy_l, grid_y_bdy_r]', [Fields_yl_List{expt}(:, 4), Fields_yr_List{expt}(:, 4)]', 'r')
        hold on
        plot(subplotax_cur, y_coord, kappa_exact(1,:), 'b:', ...
            y_coord, kappa_exact(2,:), 'r:')
        xlabel('$y$');
        run('Script_ax_curesConfig.m')

        %legend('$\kappa_{xx}$', '$\kappa_{yy}$') %, {}, {})
        %ylabel(['h = ', num2str(SamplingInterval)]);
    end
    %suptitle('Posterior Mean Fields')
    subplotax_cur = subplot(Nexpt_plot,2, 1);
    title('$u(y), v(y)$')
    run('Script_ax_curesConfig.m')
    
    subplotax_cur = subplot(Nexpt_plot,2, 2);
    title('$\kappa_{xx}(y), \kappa_{yy}(y)$')
    run('Script_ax_curesConfig.m')

    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
    % -------------------------------------------------------% 

    % Part 2: Fields at chosen points vs sampling intervals
    FigObj = figure('Name', 'vs SIv', 'NumberTitle','off', 'pos', [200 200 1200 800]);
    FigObjax_cur = ax_cures;
    legend_label = {};
    
   	% Plot pointwise diffusivity under different sampling intervals  
    % Choose 3 points to visulise
    % Assumed the point is on the boundary of the cell
    j_y000 = find(grid_y_bdy_l == 0);
    j_y025 = find(grid_y_bdy_l == 0.25);   % Ad-hoc assumed power of 100
    j_y050 = find(grid_y_bdy_l == 0.5);   % Ad-hoc assumed power of 4
    y000 = 0;
    y025 = 0.25;   % Ad-hoc assumed power of 4
    y050 = 0.5;
    U_00 = veloc_fldStruct.u(0, y000); K_00 = kappa_fldStruct.kappa(0, y000);
    U_03 = veloc_fldStruct.u(0, y025); K_03 = kappa_fldStruct.kappa(0, y025);
    U_05 = veloc_fldStruct.u(0, y050); K_05 = kappa_fldStruct.kappa(0, y050);

    % Set up storage
    FineField_array = zeros(Nexpt, 5, 3);  % 5:u,v, Kxx,Kyy,Kxy; 3: 3 points
    SIv_List = zeros(Nexpt, 1);
    for expt = 1:Nexpt
        SIv_List(expt) = DownSampling_Param_List{expt}.SamplingInterval;
        
        FineField_array(expt, :, 1) = 0.5*( Fields_yl_List{expt}(j_y000,:) + Fields_yr_List{expt}(j_y000-1,:) );
        FineField_array(expt, :, 2) = 0.5*( Fields_yl_List{expt}(j_y025,:) + Fields_yr_List{expt}(j_y025-1,:) );
        FineField_array(expt, :, 3) = 0.5*( Fields_yl_List{expt}(j_y050,:) + Fields_yr_List{expt}(j_y050-1,:) );
    end
    
    % Plot
    u_exact = [SIv_List*0+U_00(1), SIv_List*0+U_03(1), SIv_List*0+U_05(1)];
    KXX_exact = [SIv_List*0+K_00(1,1), SIv_List*0+K_03(1,1), SIv_List*0+K_05(1,1)];
    KYY_exact = [SIv_List*0+K_00(2,2), SIv_List*0+K_03(2,2), SIv_List*0+K_05(2,2)];
    
    FigObj = figure('Name', 'vs SI', 'NumberTitle','off', 'pos', [200 200 1200 800]);
    subplotax_cur = subplot(1, 3, 1);
    hplot1 = plot(subplotax_cur, SIv_List, squeeze(FineField_array(:, 1, :))');
    set(hplot1, {'Color'}, {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]});    % Specify colours of lines
    xlabel('h');
    %legend('y=0', 'y=0.25', 'y=0.5')
    hold on
    hplot2 = plot(subplotax_cur, SIv_List, u_exact, ':');
    set(hplot2, {'Color'}, {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]});
    ylim([-0.2, 1.4])
    title('$u(y)$')
    legend(hplot1, {'y=0', 'y=0.25', 'y=0.5'}, 'Location','northwest', ...
                'Orientation','vertical');
    run('Script_ax_curesConfig.m')

    subplotax_cur = subplot(1, 3, 2);
    hplot1 = plot(subplotax_cur, SIv_List, squeeze(FineField_array(:, 3, :))');
    xlabel('Sampling Interval (h)');
    hold on
    hplot2 = plot(subplotax_cur, SIv_List, KXX_exact, ':');
    legend(hplot1, {'y=0', 'y=0.25', 'y=0.5'}, 'Location','northwest', ...
                'Orientation','vertical');
    title('$K_{xx}(y)$')

    run('Script_ax_curesConfig.m')

    subplotax_cur = subplot(1, 3, 3);
    hplot1 = plot(subplotax_cur, SIv_List, squeeze(FineField_array(:, 4, :))');
    xlabel('h');
    hold on
    hplot2 = plot(subplotax_cur, SIv_List, KYY_exact, ':');
    legend(hplot1, {'y=0', 'y=0.25', 'y=0.5'}, 'Location','northwest', ...
                'Orientation','vertical');
    title('$K_{yy}(y)$')
    run('Script_ax_curesConfig.m')



 	% Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_FieldPts', '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
   	% -------------------------------------------------------% 
    
    % Part 3: Error in Diffusivity Diagnosis
    FigObj = figure('Name', 'Kxx vs SIv', 'NumberTitle','off', 'pos', [200 200 1200 800]);
    FigObjax_cur = ax_cures;
    
    FigObjax_cur = subplot(2, 1, 1)
    u_error = (abs(squeeze(FineField_array(:, 1, :)) - u_exact)./1)';
    plot(FigObjax_cur, log10(SIv_List), log10(u_error))
    xlabel('$\log(\frac{h}{1})$');
    ylabel('$\log(|\frac{\bar{u}-u}{U_0}|)$');
    xlim([-2.15, 1])
    ylim([-6, 0])
    
    FigObjax_cur = subplot(2, 1, 2)
    KXX_error = (abs(squeeze(FineField_array(:, 3, :)) - KXX_exact)./KXX_exact)';
    plot(FigObjax_cur, log10(SIv_List), log10(KXX_error))
    xlabel('$\log(\frac{h}{1})$');
    ylabel('$\log(|\frac{\bar{\kappa}_{xx}-\kappa_{xx}}{\kappa_{xx}}|)$');
    xlim([-2.15, 1])
    ylim([-4, 2.2])
    
    legend_List = {'y=0', 'y=0.25', 'y=0.5'}    
    legend(legend_List, 'Location','northwest', ...
                'Orientation','vertical');

%     FigObjax_cur = subplot(3, 1, 3)
%     KXX_error = (abs(squeeze(FineField_array(:, 4, :)) - KYY_exact)./KYY_exact)';
%     plot(FigObjax_cur, log10(SIv_List), log10(KXX_error))
%     xlabel('$\log(h)$');
%     ylabel('$\log(|\frac{\bar{\kappa}_{yy}-\kappa_{yy}}{\kappa_{yy}}|)$');
%     ylim([-4, 1.5])
    
            
    %title('log-log plot of error in $K_{XX}$');
    for i = 1:2
        FigObjax_cur = subplot(2, 1, i)
        run('Script_ax_curesConfig.m')
    end
    
    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , veloc_Profile , '_KerrorSIV', '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
    % Additional work for comparison
    DataFile_path = [ Fig_outputfolder, veloc_Profile, '_InferredFields.mat']
    save(DataFile_path, 'SIv_List', 'FineField_array', ...
            'u_exact', 'KXX_exact', 'KYY_exact', '-v7.3')
end


if strcmp(RunProfile, 'cosine_pwl_vs_pwc')
    veloc_Profile = 'cosine_pwc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    DataFile_path = [ Fig_outputfolder, veloc_Profile, '_InferredFields.mat']
    cosine_TC_pwcData = load(DataFile_path);
    
    veloc_Profile = 'cosine';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    DataFile_path = [ Fig_outputfolder, veloc_Profile, '_InferredFields.mat']
    cosine_TC_Data = load(DataFile_path);
    
    
    % Part 3: Error in Diffusivity Diagnosis
    FigObj = figure('Name', 'pwl vs pwc', 'NumberTitle','off', 'pos', [200 200 1200 800]);
    FigObjax_cur = ax_cures;
    
    DataList = [cosine_TC_Data, cosine_TC_pwcData];
    MarkerList = {'-', '-.'};
    
    for i = 1:length(DataList)
        Data = DataList(i);
        
        FigObjax_cur = subplot(2, 1, 1);
        set(FigObjax_cur, 'XDir','reverse')

        u_error = (abs(squeeze(Data.FineField_array(:, 1, :)) - u_exact)./1)';
        hplot1 = plot(FigObjax_cur, log10(SIv_List), log10(u_error), MarkerList{i});
        set(hplot1, {'Color'}, {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]});
        set(FigObjax_cur, 'XDir','reverse')

        %xlabel('$\log(\frac{h}{1})$');
        ylabel('$\log(|\frac{\bar{u}-u}{U_0}|)$');
        xlim([-2.15, 1])
        ylim([-6, 0])
        hold on
        
        FigObjax_cur = subplot(2, 1, 2);
        set(FigObjax_cur, 'XDir','reverse')

        KXX_error = (abs(squeeze(Data.FineField_array(:, 3, :)) - KXX_exact)./KXX_exact)';
        hplot2 = plot(FigObjax_cur, log10(SIv_List), log10(KXX_error), MarkerList{i});
        set(hplot2, {'Color'}, {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]});
        
        xlabel('$\log(\frac{h}{1})$');
        ylabel('$\log(|\frac{\bar{\kappa}_{xx}-\kappa_{xx}}{\kappa_{xx}}|)$');
        xlim([-2.15, 1])
        ylim([-4, 2.2])
        hold on
        
        % Estimate order of convergence for
        % large h
        small_h_filter = (log(SIv_List) < 0.5);
        SIv_List_filter = SIv_List(small_h_filter);
        
        u_error_filter = u_error(small_h_filter);        
        KXX_error_filter = KXX_error(small_h_filter);
        
        p = polyfit(log(SIv_List_filter),log(u_error_filter),1);
        disp(['u: Small h order = ', num2str(p(1))]);
        
        p = polyfit(log(SIv_List_filter),log(KXX_error_filter),1);
        disp(['Kappa: Small h order = ', num2str(p(1))]);
        
        % large h
        large_h_filter = (log(SIv_List) > 0);
        SIv_List_filter = SIv_List(large_h_filter);
        
        u_error_filter = u_error(large_h_filter);        
        KXX_error_filter = KXX_error(large_h_filter);
        
        p = polyfit(log(SIv_List_filter),log(u_error_filter),1);
        disp(['u: Large h order = ', num2str(p(1))]);
        
        p = polyfit(log(SIv_List_filter),log(KXX_error_filter),1);
        disp(['Large h order = ', num2str(p(1))]);
    end
    
    legend_List = {'y=0', 'y=0.25', 'y=0.5'};  
    legend(legend_List, 'Location','northeast', ...
                'Orientation','vertical');    
            
    %title('log-log plot of error in $K_{XX}$');
    for i = 1:2;
        FigObjax_cur = subplot(2, 1, i);
        run('Script_ax_curesConfig.m')
    end
    
    
    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, RunProfile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_KerrorSIV', '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end

if (strcmp(RunProfile, 'tg_w_mean_TC_vs_Npart') || strcmp(RunProfile, 'tg_w_mean_TC_pwc_vs_Npart'))
    %% Output histograms of posterior distributions
    Fields1theta0 = 1;
    
    Nexpt = length(theta_Stat_List);
    
    theta_mean_List = {};
    theta_vari_List = {};
    
    % Output screen
    FigObj = figure('Name','Histograms','NumberTitle','off', 'pos', [200 200 1200 800]);
    FigObjax_cur = ax_cures;
    legend_label = {};
    
    for expt = 1:Nexpt
        % Load and Process data
        MCMC_Param = MCMC_Param_List{expt};
        theta_Stat = theta_Stat_List{expt};
        DownSampling_Param = DownSampling_Param_List{expt};
        
        legend_label{end+1} = num2str(DownSampling_Param.nparticles);
        
        burnin_bgn = round(MCMC_Param.Nsteps_pc *0.25);
        %burnin_bgn = 1;

        
        % Only one cell: first 1
        theta_samples = squeeze(theta_Stat.theta_store(1, :, burnin_bgn:end))';
        [Fields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(theta_samples);
        
        
        if Fields1theta0 == 1
            theta_samples_vis = Fields(:, [1, 3, 4, 5]);
        else
            theta_samples_vis = theta_samples(:, [1, 7, 8, 9]);
        end
        
        theta_mean_List{end+1} = mean(theta_samples_vis, 1);
        theta_vari_List{end+1} = var(theta_samples_vis, 0, 1);
        
        % Plot histograms
        subplotax_cur = subplot(1,3,1);
        h11 = histogram(subplotax_cur,theta_samples_vis(:, 1), 16);
        % Edit histograms
        h11.Normalization = 'probability';
        xlabel(subplotax_cur, '$b_1$')
        %title(num2str(theta_mean_List{1}(1)))
        %xlim([0.7,0.95])
        %legend(subplotax_cur, legend_label, 'Location', 'northeast')
        hold on
        
        subplotax_cur = subplot(1,3,2);
        h2 = histogram(subplotax_cur, theta_samples_vis(:, 2), 16);
        h2.Normalization = 'probability';
        xlabel(subplotax_cur, '$\kappa_{xx}$')
        %title(num2str(theta_mean_List{1}(2)))
        %xlim([0,0.6])
        %legend(subplotax_cur, legend_label, 'Location', 'northeast')
        hold on
        
        subplotax_cur = subplot(1,3,3);
        h3 = histogram(subplotax_cur, theta_samples_vis(:, 4), 16);
        h3.Normalization = 'probability';
        xlabel(subplotax_cur, '$\kappa_{xy}$')
        %title(num2str(theta_mean_List{1}(3)))
        %xlim([0,0.4])
        legend(subplotax_cur, legend_label, 'Location', 'northeast')
        hold on
        
%         subplotax_cur = subplot(1,4,4);
%         h4 = histogram(subplotax_cur, theta_samples_vis(:, 4), 16);
%         h4.Normalization = 'probability';
%         xlabel(subplotax_cur, '$\kappa_{xy}$')
%         %title(num2str(theta_mean_List{1}(4)))
%         %xlim([0.05,0.15])
%         %legend(subplotax_cur, legend_label, 'Location', 'northeast')
%         hold on
    end
        
    % Add a vertical line for Exact Values
    %legend_label{end+1} = 'Exact';

    subplotax_cur = subplot(1,3,1);
    u_exact = veloc_fldStruct.u(0,0);
    line([u_exact(1) u_exact(1)], [0 0.3], 'Color','black', ...
        'LineStyle','--', 'LineWidth',2);
    %legend(subplotax_cur, legend_label)
    %xlim([0.75, 0.93])
    ylabel('Normalised Histogram')
    
    % Put legend onto the 4th subplot
    subplotax_cur = subplot(1,3,3);
    legend(subplotax_cur, legend_label, 'Location', 'northeast')

    
    % Exact sigma1 and phi TO BE DETERMINED

    % Configure the subplots
    for i = 1:3
        subplotax_cur = subplot(1,3,i);
        run('Script_ax_curesConfig.m')
    end

    
    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end


if strcmp(RunProfile, 'tg_w_mean_TC_vs_SIv') || strcmp(RunProfile, 'tg_w_mean_TC_pwc_vs_SIv')
    %% Output plot of inferred quantities vs Sampling Intervals
    
    Nexpt = length(theta_Stat_List);
    
    theta_mean_List = {};
    theta_vari_List = {};
    
    % Output screen
    FigObj = figure('Name', 'vs SI', 'NumberTitle','off', 'pos', [200 200 800 600]);
    FigObjax_cur = ax_cures;
    legend_label = {};
    
    for expt = 1:Nexpt
        % Load and Process data
        MCMC_Param = MCMC_Param_List{expt};
        theta_Stat = theta_Stat_List{expt};
        DownSampling_Param = DownSampling_Param_List{expt};
        
        legend_label{end+1} = num2str(DownSampling_Param.SamplingInterval);
        
        burnin_bgn = round(MCMC_Param.Nsteps_pc *0.25);
        
        % Only one cell: first 1
        theta_samples = squeeze(theta_Stat.theta_store(1, :, burnin_bgn:end));
        
        theta_mean_List{end+1} = mean(theta_samples, 2);
        theta_vari_List{end+1} = var(theta_samples, 0, 2);
    end
    
    % Post-process the mean and variance data
    mean_array = zeros(length(theta_vari_List), length(theta_vari_List{1}));
    varrel_array = zeros(length(theta_vari_List), length(theta_vari_List{1}));
    log_nsample = zeros(length(theta_vari_List), 1);
    SIv_List = zeros(length(theta_vari_List), 1);
    
    for i = 1:length(theta_vari_List)
        mean_array(i, :) = theta_mean_List{i};
        varrel_array(i, :) = theta_vari_List{i}./(theta_mean_List{i}.^2);
        
        log_nsample(i) = log(DownSampling_Param_List{i}.nparticles);
        SIv_List(i) = DownSampling_Param_List{i}.SamplingInterval;
    end
    
    [MeanFields, GradxFields, GradyFields] = convert_theta_U1K0_to_Fields(mean_array);
    
    FigObjax_cur = subplot(2,1,1)
    plot(FigObjax_cur, SIv_List, MeanFields(:, 1:2))
    legend('$b_{1}$', '$b_{2}$', ...
            'Location', 'northeast', 'orientation', 'horizontal')
    xlim([0, 64])
    ylim([0.5, 1.1])
    ylabel('Velocity')
    run('Script_ax_curesConfig.m')

    FigObjax_cur = subplot(2,1,2)
    hplot1 = plot(FigObjax_cur, SIv_List, MeanFields(:, 3:5))
    legend('$\kappa_{xx}$', '$\kappa_{yy}$', '$\kappa_{xy}$', ...
            'Location', 'northeast', 'orientation', 'horizontal')
    set(hplot1, {'Color'}, {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]});    % Specify colours of lines
    hold on
    % Add exact values
    
    
    xlim([0, 64])
    ylim([0. 0.4])
    ylabel('Diffusivity')
    xlabel('sampling interval ($h$)');
    run('Script_ax_curesConfig.m')

    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end



if strcmp(RunProfile, 'childress_soward_TC_vs_shear')
    % Compute theta_mean_List and theta_vari_List
    MCMC_NSample = MCMC_Stat_List{1}.MCMC_NSample;
    theta_mean_List = theta_Stat_List{1}.theta_sum/MCMC_NSample;
    thetasq_mean = theta_Stat_List{1}.thetasq_sum/MCMC_NSample;
    
    theta_vari_List = thetasq_mean - theta_mean_List.^2;
    
    % Loading Data
    Mesh = Mesh_Struct.Mesh;
    RectMesh_Param = Mesh_Struct.RectMesh_Param;
    
    % theta_sampling_sd = MHRW_Param_List{1}.theta_sampling_sd;
    
    NJumps_cell = MCMC_Param_List{1}.NJumps_cell;
    
    % Load TrajJumps_DA
    BTD = load(filename_BTD);
    TrajJumps_MomentGlobal = BTD.TrajJumps_DA.TrajJumps_MomentGlobal;
    JumpsCFL = BTD.TrajJumps_DA.JumpsCFL;
    
    
    grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
    NJumps_cell_k = convert_cellk_to_cellij( NJumps_cell, RectMesh_Param );
    JumpsCFL_k = convert_cellk_to_cellij( JumpsCFL, RectMesh_Param );
    
    % Crude Estimates of Velocities: by simple mean of jumps
    U_mjumps = TrajJumps_MomentGlobal.diffx(:, 1)./TrajJumps_MomentGlobal.h(:,1);
    V_mjumps = TrajJumps_MomentGlobal.diffy(:, 1)./TrajJumps_MomentGlobal.h(:,1);
    
    U_mjumps_k = convert_cellk_to_cellij( U_mjumps, RectMesh_Param );
    V_mjumps_k = convert_cellk_to_cellij( V_mjumps, RectMesh_Param );
    
    % Crude Estimates of Diffusivity: by simple mean of jumps
    Kxx_jumps = TrajJumps_MomentGlobal.diffx(:, 2)./(2*TrajJumps_MomentGlobal.h(:,1));
    Kyy_jumps = TrajJumps_MomentGlobal.diffy(:, 2)./(2*TrajJumps_MomentGlobal.h(:,1));
    
    Kxx_jumps_k = convert_cellk_to_cellij(Kxx_jumps, RectMesh_Param );
    Kyy_jumps_k = convert_cellk_to_cellij(Kyy_jumps, RectMesh_Param );
    
    
    
    % Output screen
    % Configure figure (Currently not in use)
    fig_xI = 0.05; fig_yI = 0.05;
    fig_dx = 0.9; fig_dy = 0.9;

    % Range of plot
    U_ub = 1.5*max_cur([abs(U_mjumps(:)); abs(V_mjumps(:))]);
    U_range = [-U_ub U_ub];
    U_CLV = linspace(U_range(1), U_range(2), 20);
    
   
    % Obtain the 2D Field Arrays
    [MeanFields_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_mean_List, Mesh_Struct);
    
    % 1. Evaluate Exact Fields
    u_exact = zeros(size(grid_vis_f.X));
    v_exact = zeros(size(grid_vis_f.X));
    sigma1_exact = zeros(size(grid_vis_f.X));
    sigma2_exact = zeros(size(grid_vis_f.X));
    phi_exact = zeros(size(grid_vis_f.X));
    
    for i = 1:size(grid_vis_f.X)
        for j = 1:size(grid_vis_f.X)
            x_ij = grid_vis_f.X(i,j);
            y_ij = grid_vis_f.Y(i,j);
            
            U_ij = veloc_fldStruct.u(x_ij, y_ij);
            K_ij = kappa_fldStruct.kappa(x_ij, y_ij);
            
            u_exact(i,j) = U_ij(1);
            v_exact(i,j) = U_ij(2);
            
            [K_sigma1, K_sigma2, K_phi]  = Kcart_to_Kpolar(K_ij(1,1), K_ij(2,2), K_ij(1,2));
            sigma1_exact(i,j) = K_sigma1;
            sigma2_exact(i,j) = K_sigma2;
            phi_exact(i,j) = K_phi;
        end
    end
    
    % 2a. Crude Estimate of Velocity Fields
    FineGrid_Resolution = 100;
    ZeroGrad = zeros(size(U_mjumps_k));
    [U_mjumps_f, X_f, Y_f] = griddata_interpolation_pwl(U_mjumps_k, ZeroGrad, ZeroGrad, RectMesh_Param, FineGrid_Resolution);
    [V_mjumps_f, X_f, Y_f] = griddata_interpolation_pwl(V_mjumps_k, ZeroGrad, ZeroGrad, RectMesh_Param, FineGrid_Resolution);
    
    FigObj = figure('Name','Velocity','NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    contourf_vis_data = {u_exact, v_exact; ...
        U_mjumps_f, V_mjumps_f; ...
        MeanFields_f.u, MeanFields_f.v; }';
    ContourLineValue_data = {U_CLV, U_CLV; U_CLV, U_CLV; U_CLV, U_CLV}';
    cax_curis_data = {U_range, U_range; U_range, U_range; U_range, U_range}';
    colormap_data = {'jet', 'jet'; 'jet', 'jet'; 'jet', 'jet'}';
    
    title_data = {'u: Exact', 'v: Exact'; 'w/o shear', 'w/o shear'; 'w shear', 'w shear'}';
    subplot_param = struct('Position', [fig_xI fig_yI fig_dx fig_dy], ...
        'ContourLineValue_data', ContourLineValue_data, ...
        'cax_curis_data', cax_curis_data, 'title_data', title_data, ...
        'colormap_data', colormap_data);
    
    contourf_subplot(grid_vis_f, contourf_vis_data, subplot_param);
    
    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_veloc.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
    
    % 2b. Crude Estimate of Diffusivity in Cartesian Coordinates
    [Kxx_exact, Kyy_exact, Kxy_exact]  = Kpolar_to_Kcart_vectorised(sigma1_exact, sigma2_exact, phi_exact);
    [Kxx_mean, Kyy_mean, Kxy_mean]  = Kpolar_to_Kcart_vectorised(MeanFields_f.K_sigma1, MeanFields_f.K_sigma2, MeanFields_f.K_phi);
    Kxx_range = [min(Kxx_exact(:)), max_cur(Kxx_exact(:))];
    Kyy_range = [min(Kyy_exact(:)), max_cur(Kyy_exact(:))];
    Kxy_range = [min(Kxy_exact(:)), max_cur(Kxy_exact(:))];
    Kxx_CLV = linspace(Kxx_range(1), Kxx_range(2), 20);
    Kyy_CLV = linspace(Kyy_range(1), Kyy_range(2), 20);
    Kxy_CLV = linspace(Kxy_range(1), Kxy_range(2), 20);
    
    [Kxx_f, X_f, Y_f] = griddata_interpolation_pwl(Kxx_jumps_k, Kxx_jumps_k*0, Kxx_jumps_k*0, RectMesh_Param, 100);
    [Kyy_f, X_f, Y_f] = griddata_interpolation_pwl(Kyy_jumps_k, Kyy_jumps_k*0, Kyy_jumps_k*0, RectMesh_Param, 100);
    
    FigObj = figure('Name','Exact Diffusivity (Cart, Publish)','NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    contourf_vis_data = {Kxx_exact, Kyy_exact; ...
        Kxx_f, Kyy_f; ...
        Kxx_mean, Kyy_mean;}';
    ContourLineValue_data = {Kyy_CLV, Kyy_CLV; Kyy_CLV, Kyy_CLV; Kyy_CLV, Kyy_CLV}';
    cax_curis_data = {Kyy_range, Kyy_range; Kyy_range, Kyy_range; Kyy_range, Kyy_range}';
    colormap_data = {'jet', 'jet'; 'jet', 'jet'; 'jet', 'jet'}';
    
    title_data = {'K_{xx}: Exact', 'K_{yy}: Exact'; 'w/o shear', 'w/o shear'; 'w/ shear', 'w/ shear'}';
    
    subplot_param = struct('Position', [fig_xI fig_yI fig_dx fig_dy], ...
        'ContourLineValue_data', ContourLineValue_data, ...
        'cax_curis_data', cax_curis_data, 'title_data', title_data, ...
        'colormap_data', colormap_data);
    
    contourf_subplot(grid_vis_f, contourf_vis_data, subplot_param);
    
    % Save the figure
    ImageFormat = 'epsc';
    Fig_outputfolder = [DataFolder_Publication, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , '_kappa.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
    
end