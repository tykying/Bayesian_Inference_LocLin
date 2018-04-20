%close all
addpath('../Utilities/')
addpath('../StructConversion/')
addpath('../StatCompute/')
addpath('../Scripts/')
addpath('./PlotUtilities/')
addpath('../PlotUtilities/')

close all
%load('/data/tying/BayesianInference/BayesSampleStat/Layer_2/QGM2_DStemp/BayesSampleStat_QGM2_DStemp_DSr1_h128_Idx16_LocLinPolar_Nsteps1048576_nTrial0.mat');
%veloc_Profile = 'QGM2_DStemp';
%RunProfile = 'QGM2_DStemp';
%nTrial = 0;

% Want to load specific file?
% Parameters to control
% veloc_Profile = 'linear_shear';
% DS_rate = 1;
% SamplingInterval_vis = 1;
% Nx_cell_ARG = 1;
% Nsteps_pc = 4096;
% InfScheme = 'Local';
% Layer = 2;
% nTrial = 1;
% DataTypeString = 'BayesSampleStat';
% SamplingInterval_vis_List = [2 4 8 16 32 64 128 256];
%
% NSIv = length(SamplingInterval_vis_List);
% NTrial = 10;
%
% theta_mean_ijk = zeros(NSIv, NTrial, 9);
% error_ij = zeros(NSIv, NTrial);
%
%
% NTrial = 1;
% SamplingInterval_vis_List = [4];
% SamplingInterval_vis = 16;


%run('Script_Filenames');
%load(filename_BSS);
theta_Stat = theta_Stat_List{nTrial};

%% Visualise all particle positions
plotTraj = 0;
if plotTraj
    FigObj = figure('Name','Data','NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    subplot(2,2,1)
    x_vis = TrajData.x;
    y_vis = TrajData.y;
    plot(x_vis, y_vis)
    title('Trajectories')
    
    subplot(2,2,2)
    scatter(x_vis(end,:), y_vis(end,:))
    title('Final Positions of Particles')
    
    subplot(2,2,3)
    scatter(x_vis(1,:), y_vis(1,:))
    title('Initial Positions of Particles')
    
    subplot(2,2,4)
    h1 = histogram(x_vis(end,:));
    h1.Normalization = 'probability';
    hold on
    h2 = histogram(y_vis(end,:));
    h2.Normalization = 'probability';
    title('Distributions of Final Positions of Particles')
end


%% Check MC stationary
nTrial_str = num2str(nTrial);

Ntheta_store = MCMC_Stat.Ntheta_store;
Ntheta_store_ZIP = size(theta_Stat.theta_store, 3);
SamplingInterval_vis = DownSampling_Param.SamplingInterval_vis;
InfScheme = MCMC_Param.InfScheme;
Nx_cell_ARG = Mesh_Struct.RectMesh_Param.Nx_cell;

% Conversion from theta to theta_Fields (Cartesian)
% theta_Fields = theta_Stat.theta_store(:,:,round(Ntheta_store_ZIP*0.5):end);
[Ncells, NVars, MCMC_NSample] = size(theta_Stat.theta_store);

theta_Fields = convert_theta_to_theta_Fields(theta_Stat.theta_store);


%% Convert theta into theta_Fields
% Only check convergence of 1st and 2nd components
cmp_List = [1,3]; 
Ncmp_TF = length(cmp_List);

MCStat_TF_cmp = zeros(Ncells, Ncmp_TF);
fprintf('Checking Markov Chain Stationarity. \n');

for cell_k = 1:Ncells
    for cmp_ind = 1:Ncmp_TF
        cmp = cmp_List(cmp_ind);
        
        MC_Data = theta_Fields(cell_k, cmp, :);
        MCStat_TF_cmp(cell_k, cmp_ind) = check_MarkovChain_stationarity(MC_Data);
    end
end

% Locate the cells that diverged
cell_k_List_OneDiverged = find(sum(MCStat_TF_cmp, 2) == 1);
cell_k_List_TwoDiverged = find(sum(MCStat_TF_cmp, 2) == 0);

cell_k_List_Diverged = find(sum(MCStat_TF_cmp, 2) < Ncmp_TF);

%%
tic;
NSessions = 4;
if NSessions > 1
    theta_Fields_Ses = cell(1,NSessions);
    for Ses = 1:NSessions
        Ses_bgn = round(MCMC_NSample*(Ses-1)/NSessions)+1;
        Ses_end = round(MCMC_NSample*Ses/NSessions);
        fprintf('Ses: %d, in samples (%d, %d). \n', Ses, Ses_bgn, Ses_end);
        
        theta_Fields_Ses{Ses} = theta_Fields(:, :, Ses_bgn:Ses_end);
    end    
end
DT_conKpolar = toc;
fprintf('Finished session MC in %f sec. \n', DT_conKpolar);


%% MAP Estimates
theta_MAP = zeros(Ncells, NVars);
% burnin = round(MCMC_NSample*0.5);

tic
logPost_MAP_est = 1;
if logPost_MAP_est == 1       
    % Read from logPost
    theta_MAP = convert_theta_to_theta_Fields(theta_Stat.theta_MAP);

else
    % Read distribution
    SamplesStored_ARRAY = theta_Fields;

    parfor cell_k = 1:Ncells
        if mod(log2(cell_k), 4) == 0
            fprintf('Working on %d. \n', cell_k);
        end
        
        theta_MAP_cellk = zeros(1, NVars);
        SamplesStored_cellk = squeeze(SamplesStored_ARRAY(cell_k,:,:));
        
        for cmp = 1:NVars
            SamplesStored = SamplesStored_cellk(cmp, :);
            theta_MAP_cellk(cmp) = estimate_mode(SamplesStored);
        end
        
        theta_MAP(cell_k, :) = theta_MAP_cellk;
    end
end


Dt_EstimateMAP = toc;
fprintf('Finsihed Sesing MC in %f sec. \n', Dt_EstimateMAP);

%% Plot the Markov Chain over iterations
plotMC = 0;
if plotMC == 1
    assert(MCMC_NSample == (size(Fields, 3)));
    chain_len = size(Fields, 3);
    chain = 1:chain_len';

    % Chain of velocity
    FigObj = figure('Name', [nTrial_str, ': Mean Velocity'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ind_List = [1:2];
    title_List = {'u', 'v'};
    for i = 1:length(title_List)
        subplot(2,1,i)
        data = squeeze(theta_Fields(1:4:end, ind_List(i), :));
        plot(chain, data')
        title(title_List{i})
        xlim([1, chain_len]);
    end
    
    
    % Chain of velocity gradient
    FigObj = figure('Name', [nTrial_str, ': Velocity Gradients'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ind_List = [3:5];
    title_List = {'A_{11}', 'A_{12}', 'A_{21}'};
    
    for i = 1:length(title_List)
        subplot(3,1,i)
        data = squeeze(theta_Fields(1:4:end, ind_List(i), :));
        plot(chain, data')
        title(title_List{i})
        xlim([1, chain_len]);
    end
    
    % Chain of diffusivity
    FigObj = figure('Name', [nTrial_str, ': Diffusivity'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ind_List = [7:9];
    title_List = {'sigma_1', 'sigma_2', 'phi'};
    Fldind_List = [3:5];
    Fldtitle_List = {'K_{xx}', 'K_{yy}', 'K_{xy}'};
    
    for i = 1:length(title_List)
        subplot(3,2,2*i-1)
        data = squeeze(theta_Fields(1:4:end, ind_List(i), :));
        plot(chain, data')
        title(title_List{i})
        xlim([1, chain_len]);

        subplot(3,2,2*i)
        data = squeeze(Fields(1:4:end, Fldind_List(i), :));
        plot(chain, data')
        title(Fldtitle_List{i})
        xlim([1, chain_len]);
    end
end

%% Visualise MCMC Histogram
cellij_to_k = Mesh_Struct.RectMesh_Param.cellij_to_k;
cellk_to_ij = Mesh_Struct.RectMesh_Param.cellk_to_ij;

plotLowAccRat = 1;
plotLowAsp= 0;
plotLowSigma2 = 0;
plotCellsSelected= 1;
plotManuel = 0;
plotWeirdAccRat = 0;
plotOneDiverged = 0;
plotTwoDiverged = 0;
if plotLowAccRat == 0
    % Form cell_k_vis
    cell_k_AccRatlow = find(MCMC_Stat.acceptance_ratio(:, 1) < 0.05);
    cell_k_vis = cell_k_AccRatlow';
end
if plotWeirdAccRat == 1
    % Form cell_k_vis
    cell_k_WeirdAccRat = find((MCMC_Stat.acceptance_ratio(:, 1)-0.23) > 0.15);
    cell_k_vis = cell_k_WeirdAccRat';
end
if plotLowAsp == 1
    % Form cell_k_vis
    cell_k_LowAR = find(theta_MAP(:, 7)./theta_MAP(:, 8) > 3);
    cell_k_vis = cell_k_LowAR';
end
if plotCellsSelected == 1
    % Form cell_k_vis
    cell_i_List = [3, 9, 13];
    cell_j_List = [3, 9, 13];

    [cell_i_Selected, cell_j_Selected] = meshgrid(cell_i_List, cell_j_List);
    
    cell_k_vis = zeros(1, length(cell_i_Selected(:)) );
    for cell_ind = 1:length(cell_k_vis)
        cell_k_vis(cell_ind) = cellij_to_k(cell_i_Selected(cell_ind), cell_j_Selected(cell_ind));
    end
end
if plotLowSigma2 == 1
    % Form cell_k_vis
    cell_k_LowSigma2 = find(theta_MAP(:, 8) < 4);
    cell_k_vis = cell_k_LowSigma2';
end
if plotOneDiverged == 1
    cell_k_vis = cell_k_List_OneDiverged';
end
if plotTwoDiverged == 1
    cell_k_vis = cell_k_List_TwoDiverged';
end
if plotManuel == 1
    cell_k_vis = [cellij_to_k(8, 10), cellij_to_k(4, 10)];
end


% Visualisation
close all;
NCells_vis = length(cell_k_vis);

Title_List = {'u', 'v', 'd_x u', 'd_y u', 'd_x v', 'd_y v', '\sigma_1', '\sigma_2', '\phi'};


FigObj = {};
for cell_ind = 1:NCells_vis
    cell_k = cell_k_vis(cell_ind);
    cellij = cellk_to_ij(cell_k_vis(cell_ind));
    
    FigObj{cell_ind} = figure('Name', sprintf('i: %d, j: %d; k: %d', cellij(1), cellij(2), cell_k) ,'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    histfig = cell(NCells_vis, NVars, NSessions+1);
    for cmp = 1:NVars
        ax = subplot(3,3, cmp);
        cell_k = cellij_to_k(cellij(1), cellij(2));
        
        hold(ax, 'on')
        if NSessions > 1
            for Ses = 1:NSessions
                histfig{cell_ind, cmp, Ses} = histogram(ax, theta_Fields_Ses{Ses}(cell_k, cmp, :), 64);
                histfig{cell_ind, cmp, Ses}.Normalization = 'pdf';
            end
        end
        
        histfig{cell_ind, cmp, NSessions+1} = histogram(ax, theta_Fields(cell_k, cmp, :), 64);
        histfig{cell_ind, cmp, NSessions+1}.Normalization = 'pdf';
        title(Title_List{cmp})
    end
end

%% Visualise Fields
% Loading Data
Mesh = Mesh_Struct.Mesh;
RectMesh_Param = Mesh_Struct.RectMesh_Param;
cellij_to_k = Mesh_Struct.RectMesh_Param.cellij_to_k;
cellk_to_ij = Mesh_Struct.RectMesh_Param.cellk_to_ij;

theta_SSD = MHRW_Param.theta_SSD;

NJumps_cell = MCMC_Param.NJumps_cell;

grid_vis = setup_grid_vis(Mesh, RectMesh_Param);

% Choose whether or not to plot
plotJumpsData = 0;
plottheta_ssd_Init = 1;
plotFieldsMoments = 0;
plotMAPFields = 1;
plotMCStat_TF = 1;
plotMAPFields_f = 1;

if plotJumpsData == 1
    load(filename_BTD);
    
    TrajJumps_MomentGlobal = TrajJumps_DA.TrajJumps_MomentGlobal;
    JumpsCFL = TrajJumps_DA.JumpsCFL;
    
    % Configure figure
    fig_xI = 0.05; fig_yI = 0.05;
    fig_dx = 0.9; fig_dy = 0.9;
    
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
    
    X_Ind = 1:length(grid_vis.X(:,1));
    Y_Ind = 1:length(grid_vis.Y(1,:));
    
    FigObj = figure('Name','Basic','NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    subplot(2, 2, 1)
    %plot(theta_Stat.logPost_ts(:,:, 1:1024:end));
    title('Log-Post')
    subplot(2, 2, 2)
    imagesc(X_Ind, Y_Ind, NJumps_cell_k');
    axis('xy')
    %pcolor(grid_vis.X, grid_vis.Y, log10(NJumps_cell_k))
    title('log10-NJumps')
    colorbar
    subplot(2, 2, 3)
    imagesc(X_Ind, Y_Ind, JumpsCFL_k(:,:,1)');
    axis('xy')
    %pcolor(grid_vis.X, grid_vis.Y, JumpsCFL_k(:,:,1))
    title('JumpsCFL: x')
    colorbar
    subplot(2, 2, 4)
    imagesc(X_Ind, Y_Ind, JumpsCFL_k(:,:,2)');
    axis('xy')
    %pcolor(grid_vis.X, grid_vis.Y, JumpsCFL_k(:,:,2))
    title('JumpsCFL: y')
    colorbar
end


if plotMCStat_TF == 1
    % Sampler Parameters
%     MCStat_TF_ij = convert_cellk_to_cellij(MCStat_TF_cmp, RectMesh_Param);
%     FigObj = figure('Name',[nTrial_str, ': MCStat_TF_ij'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
%     ax_sp = plot_dataij_k(grid_vis, MCStat_TF_ij);
    
    % Mask the divergent cells
    theta_MAP_MASK = theta_MAP;
    theta_MAP_MASK(cell_k_List_Diverged, :) = NaN;
    %cell_k_List_Diverged_old = cell_k_List_Diverged;

    Fields_MAP_MASK_CP = convert_cellk_to_cellij(theta_MAP_MASK, RectMesh_Param);
    FigObj = figure('Name',[nTrial_str, ': theta_MAP_MASK'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ax_sp = plot_dataij_k(grid_vis, Fields_MAP_MASK_CP);
    
    %%
    ImageFormat = 'png';
    DataFolder_Figure = '/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Figure_Generated/';
    
    Fig_outputfolder = [DataFolder_Figure, veloc_Profile, '/'];
    [status, msg, msgID] = mkdir(Fig_outputfolder);
    
    Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis), '_', InfScheme, '_', num2str(nTrial)];
    
    Fig_outputpath = [ Fig_outputfolder , RunProfile , 'theta_MAP_MASKED', Figure_param, '.', ImageFormat(1:3)];
    saveas(FigObj, Fig_outputpath, ImageFormat);
end

if plottheta_ssd_Init == 1
    % Sampler Parameters
    theta_Ssd = convert_cellk_to_cellij(theta_SSD, RectMesh_Param);
    FigObj = figure('Name',[nTrial_str, ': Theta Ssd'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ax_sp = plot_dataij_k(grid_vis, theta_Ssd);
    
    Fields_InitCP = convert_cellk_to_cellij(squeeze(theta_Fields(:,:,1)), RectMesh_Param);
    FigObj = figure('Name',[nTrial_str, ': Init theta_Fields_Ses'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ax_sp = plot_dataij_k(grid_vis, Fields_InitCP);
    
    Fields_AccRat = convert_cellk_to_cellij(abs(MCMC_Stat.acceptance_ratio-0.23), RectMesh_Param);
    Fields_AccRat = convert_cellk_to_cellij(MCMC_Stat.acceptance_ratio, RectMesh_Param);
    %Fields_AccRat(Fields_AccRat < 0.15) = 0;
    FigObj = figure('Name',[nTrial_str, ': Acceptance Ratio'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ax_sp = plot_dataij_k(grid_vis, Fields_AccRat);
end

if plotFieldsMoments == 1
    %% Moments of posterior
    theta_Fields_Moments = zeros(Ncells, NVars, 4);
    
    % Summary
    theta_Fields_Moments(:, :, 1) = mean(theta_Fields, 3);
    theta_Fields_Moments(:, :, 2) = var(theta_Fields, 0, 3);
    theta_Fields_Moments(:, :, 3) = skewness(theta_Fields, 0, 3);
    theta_Fields_Moments(:, :, 4) = kurtosis(theta_Fields, 0, 3);
    
    for moment = 1:4
        Fields_Moments = convert_cellk_to_cellij(theta_Fields_Moments(:,:,moment), RectMesh_Param);
        FigObj = figure('Name',[nTrial_str, ': ', num2str(moment), ' theta_Fields_Moments'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        ax_sp = plot_dataij_k(grid_vis, Fields_Moments);
    end
end

if plotMAPFields == 1
    Fields_MAPCP = convert_cellk_to_cellij(theta_MAP, RectMesh_Param);
    FigObj = figure('Name',[nTrial_str, ': theta_MAP'],'NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
    ax_sp = plot_dataij_k(grid_vis, Fields_MAPCP);
end
%%
if plotMAPFields_f == 1
    [theta_MAP_ij_f, grid_vis_f] = convert_theta_Fields_into_FineFields(theta_MAP, Mesh_Struct);

    Fig_outputpath_lbl = {'_Uf', '_Vf'};
    for lbl = 1:2
        FigObj = figure('Name', 'Fields_UV', 'NumberTitle','on', 'pos', [200 200 800 800]);
        ax_cur = gca;
        
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
        
        ImageFormat = 'png';
        DataFolder_Figure = '/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Figure_Generated/';

        Fig_outputfolder = [DataFolder_Figure, veloc_Profile, '/'];
        [status, msg, msgID] = mkdir(Fig_outputfolder);
        
        Figure_param = ['_Idx', num2str(Nx_cell_ARG), '_h', num2str(SamplingInterval_vis), '_', InfScheme, '_', num2str(nTrial)];
        
        Fig_outputpath = [ Fig_outputfolder , RunProfile , Fig_outputpath_lbl{lbl}, Figure_param, '.', ImageFormat(1:3)];
        saveas(FigObj, Fig_outputpath, ImageFormat);
    end
end