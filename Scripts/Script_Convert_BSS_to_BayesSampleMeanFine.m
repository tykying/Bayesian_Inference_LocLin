%clear all; close all;

addpath('../Utilities/')
addpath('../StructConversion/')
addpath('../StatCompute/')
addpath('../Scripts/')


% Parameters in sliders (Integer-valued, from 1 to max)
Layer_List = [2];
Nx_cell_ARG_List = [32,];
SamplingInterval_vis_List = [1, 2, 4, 8, 16, 32, 48, 64, 80, 96, 112, 128];
%SamplingInterval_vis_List = SIv_in;
DS_rate_List = [1];
Nsteps_pc_List = [4096];
nTrial_List = [1:2];

% Parameters in popup menu
InfScheme_List = {'Local'};
veloc_Profile_List = {'QGM2_DSpart'};

for i_Nx = 1:length(Nx_cell_ARG_List)
for j_hv = 1:length(SamplingInterval_vis_List)
for m_vp = 1:length(veloc_Profile_List)
for n_is = 1:length(InfScheme_List)
for p_NS = 1:length(Nsteps_pc_List)
for q_LI = 1:length(Layer_List)
for r_Dr = 1:length(DS_rate_List)
for s_nT = 1:length(nTrial_List)
% Read the arguments
Nx_cell_ARG = Nx_cell_ARG_List(i_Nx);
SamplingInterval_vis = SamplingInterval_vis_List(j_hv);
veloc_Profile = veloc_Profile_List{m_vp};
InfScheme = InfScheme_List{n_is};
Nsteps_pc = Nsteps_pc_List(p_NS);
Layer = Layer_List(q_LI);
DS_rate = DS_rate_List(r_Dr);
nTrial = nTrial_List(s_nT);

% Display Current State:
str = sprintf('Nx_cell= %d, h= %d, %s, %s, Nsteps_pc= %d, Layer= %d, DS_rate= %d, nTrial= %d', ...
    Nx_cell_ARG, SamplingInterval_vis, veloc_Profile, InfScheme, ...
    Nsteps_pc, Layer, DS_rate, nTrial);
disp(str)


%% Conversion from BSS to MSMF

% Check whether the target exists, not the source exists
DataTypeString = 'BayesSampleMeanFine';
run('Script_Filenames');

if (exist(filename_BSMF, 'file') == 2)
    disp('File Already converted');
else
    disp('Conversion under progress.');
    % Load the corresponding file
    DataTypeString = 'BayesSampleStat';
    run('Script_Filenames');
    
    BayesSampleStat = load(filename_BSS);
    
    DownSampling_Param = BayesSampleStat.DownSampling_Param;
    
    % Compute theta_mean and theta_vari
    % N.B. Combining sets of data with same mean and variance, even with
    % different sample size, will not change the combined same mean and var
    MCMC_NSample = BayesSampleStat.MCMC_Stat.MCMC_NSample;
    theta_mean = BayesSampleStat.theta_Stat.theta_sum/MCMC_NSample;
    thetasq_mean = BayesSampleStat.theta_Stat.thetasq_sum/MCMC_NSample;
    
    % Evaluate theta_vari
    theta_vari = thetasq_mean - theta_mean.^2;
    
    % Evaluate theta_Final and theta_Init
    theta_store = BayesSampleStat.theta_Stat.theta_store;
    % Check if samples of theta has been stored
    if length(size(BayesSampleStat.theta_Stat.theta_store))==3
        theta_Final = squeeze(theta_store(:,:, end));
    else
        theta_Final = theta_store;
    end
    
    theta_Init = BayesSampleStat.theta_Init;
    
    % Posterior time-series
    logPost_ts = BayesSampleStat.theta_Stat.logPost_ts;
    
    
    % Obtain the 2D Field Arrays
    Mesh_Struct = BayesSampleStat.Mesh_Struct;
    
    [MeanFields_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_mean, Mesh_Struct);
    [VariFields_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_vari, Mesh_Struct);
    
    [theta_Init_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_Init, Mesh_Struct);
    [theta_Final_f, grid_vis_f] = convert_theta_U1K0_into_FineFields(theta_Final, Mesh_Struct);
    
    % [Ad-hoc]: Assume Not Multiple-valued Field
    % i.e. fine grid not lying on the mesh boundary of the coarse grid
    assert(size(MeanFields_f.u, 3) == 1);  
    
    
    %% Output everything to BayesSampeMeanFine 
    DataTypeString = 'BayesSampleMeanFine';
    run('Script_Filenames');
    [status,msg,msgID] = mkdir(filepath);
    save(filename_BSMF, 'grid_vis_f', 'MeanFields_f', 'VariFields_f', ...
        'theta_Init_f', 'theta_Final_f', ...
        'MCMC_NSample', 'logPost_ts', 'Mesh_Struct', 'DownSampling_Param', ...
        'Nx_cell_ARG', 'SamplingInterval_vis', 'veloc_Profile', ...
        'InfScheme', 'Nsteps_pc', 'Layer', 'DS_rate');
end


%% Load and post-process Transition Density
DataTypeString = 'BinnedTrajTransMat';
run('Script_Filenames');

% BTTM: from diffusivity diagnosis
if (exist(filename_BTTM, 'file') == 2)
    disp(['File does exist: ', filename_BTTM])
    BinnedTrajTransMat = load(filename_BTTM);
    
    if isfield(BinnedTrajTransMat,'grid_vis_f') == 1
        disp(['File Already Appended: ', filename_BTTM])
    else
        disp(['Now Appending the Transition Density ', filename_BTTM])
        % Unfold the variables
        Mesh_Struct = BinnedTrajTransMat.Mesh_Struct;
        RectMesh_Param = BinnedTrajTransMat.Mesh_Struct.RectMesh_Param;
        DownSampling_Param = BinnedTrajTransMat.DownSampling_Param;
        traj_fullpath = BinnedTrajTransMat.traj_fullpath;
        
        TransMat = BinnedTrajTransMat.TransMat;  % To process
        
        % Sum along row i = total number of Jumps starting from cell_i
        TrajJumps_celli = sum(TransMat, 2);
        
        % Reduce it into Prob Matrix: Pij = Pr(j|i); sum along each row=1
        TransDen = TransMat./TrajJumps_celli;
        
        % Power method <=> Physical method to obtain the stationary distribution (pretty stupid...)
        % pi*Pij = pi
        pi_init = rand(1, size(TransMat, 1));
        pi_old = pi_init;
        for iteration = 1:50000
            pi = pi_old*TransDen;
            
            %disp(max(pi./pi_old))
            pi_old = pi;
        end
        
        pi_Stat = pi/sum(pi);
        
        % Proportion of particles remained in the same cell = diagonal terms of Pij
        RatioRemain = diag(TransDen);
        
        % Proportion of particles Transit to neighbouring cells (5-points
        % stencil)
        RatioRemainNeigh = zeros(size(RatioRemain));
        for cell_k = 1:length(RatioRemainNeigh)
            % Brainless method: first convert back to 2D array and only care +-1
            % along i and j directions
            
            % Note: basically it is just the sum of some terms along the row of the
            % sparse
            TransDen_from_k = TransDen(cell_k, :);
            %TransDen_from_k_xy = convert_cellk_to_cellij( TransDen_from_k, RectMesh_Param );
            
            % Note:
            % TransDen_from_k_xy(cell_x, cell_y) store the prob from cell_k to (cell_x, cell_j)
            
            cell_xy = RectMesh_Param.cellk_to_ij(cell_k);
            cell_x = cell_xy(1);
            cell_y = cell_xy(2);
            
            % Ensure not out of boundaries; 5-points stencil
            if (cell_x+1) <= RectMesh_Param.Nx_cell
                cell_l = RectMesh_Param.cellij_to_k(cell_x+1, cell_y);
                RatioRemainNeigh(cell_k) = RatioRemainNeigh(cell_k) + TransDen_from_k(cell_l);
            end
            
            if (cell_x-1) >= 1
                cell_l = RectMesh_Param.cellij_to_k(cell_x-1, cell_y);
                RatioRemainNeigh(cell_k) = RatioRemainNeigh(cell_k) + TransDen_from_k(cell_l);
            end
            
            if (cell_y+1) <= RectMesh_Param.Ny_cell
                cell_l = RectMesh_Param.cellij_to_k(cell_x, cell_y+1);
                RatioRemainNeigh(cell_k) = RatioRemainNeigh(cell_k) + TransDen_from_k(cell_l);
            end
            
            if (cell_y-1) >= 1
                cell_l = RectMesh_Param.cellij_to_k(cell_x, cell_y-1);
                RatioRemainNeigh(cell_k) = RatioRemainNeigh(cell_k) + TransDen_from_k(cell_l);
            end
        end
        
        
        % Conversion to 2D array form
        [ pi_Stat_ij ] = convert_cellk_to_cellij( pi_Stat', RectMesh_Param );
        [ RatioRemain_ij] = convert_cellk_to_cellij( RatioRemain, RectMesh_Param );
        [ RatioRemainNeigh_ij] = convert_cellk_to_cellij( RatioRemainNeigh, RectMesh_Param );
        
        
        FineGrid_Resolution = 100;
        [pi_Stat_f, X_f, Y_f] = griddata_interpolation_pwl(pi_Stat_ij, zeros(size(pi_Stat_ij)), zeros(size(pi_Stat_ij)), RectMesh_Param, FineGrid_Resolution);
        [RatioRemain_f, X_f, Y_f] = griddata_interpolation_pwl(RatioRemain_ij, zeros(size(RatioRemain_ij)), zeros(size(RatioRemain_ij)), RectMesh_Param, FineGrid_Resolution);
        [RatioRemainNeigh_f, X_f, Y_f] = griddata_interpolation_pwl(RatioRemainNeigh_ij, zeros(size(RatioRemainNeigh_ij)), zeros(size(RatioRemainNeigh_ij)), RectMesh_Param, FineGrid_Resolution);
        
        % Set up grid_vis_f; something is missing
        [Mesh_f, RectMesh_Param_f] = setup_RectMeshStruct(size(X_f, 1), size(Y_f, 2), RectMesh_Param.SpDis_ranges_min, RectMesh_Param.SpDis_ranges_max);
        
        grid_vis_f = setup_grid_vis(Mesh_f, RectMesh_Param_f);
        
        % Output: by appending the BTTM file
        [status,msg,msgID] = mkdir(filepath);
        save(filename_BTTM, 'TransMat', 'grid_vis_f', ...
            'pi_Stat_f', 'RatioRemain_f', 'RatioRemainNeigh_f', ...
            'pi_Stat_ij', 'RatioRemain_ij', 'RatioRemainNeigh_ij', ...
            'DownSampling_Param', 'Mesh_Struct', 'traj_fullpath');
    end
else
    disp(['BTTM File DOES NOT EXIST: ', filename_BTTM])
end

end
end
end
end
end
end
end
end