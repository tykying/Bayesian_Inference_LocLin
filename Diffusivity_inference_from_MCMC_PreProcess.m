%% Request paths
addpath('./Utilities/')
addpath('./StructConversion/')
addpath('./StatCompute/')
addpath('./MCMC/')

%% Require variables:
% 1) Layer; 2) veloc_Profile
% 3) SamplingInterval_vis; 4) DS_rate; 5) Nx_cell_ARG
% 6) InfScheme_iterator; 7) Nsteps_pc
% 8a) u_scale_i; 8b) u_angle;  [For TG]

disp(datetime('now'))
disp(['veloc_Profile = ', veloc_Profile, '; Nx_cell = ', num2str(Nx_cell_ARG), '; SamplingInterval_vis = ', num2str(SamplingInterval_vis)]);
disp(['DS_rate = ', num2str(DS_rate), '; InfScheme = ', InfScheme, '; Nsteps_pcv = ', num2str(Nsteps_pcv)]);

%% Loading And Initialise Data
run('./Scripts/Script_LoadTrajData');
disp('Finished Loading Raw Trajectory Data.')

% Info about the trajectories
[Nts_RAW, nparticles_RAW] = size(x);
NJumps_RAW = (Nts_RAW-1)*nparticles_RAW;

disp('Full Trajectories: ')
disp(['NJumps_RAW = ', num2str(NJumps_RAW), '; Nts_RAW = ', num2str(Nts_RAW), '; nparticles_RAW = ', num2str(nparticles_RAW)]);

%% DownSampling Trajectories
% Sub-sampling Trajectories
run('./Scripts/Script_SubSampling_Traj');
disp('Finished SubSampling Trajectories.')
%clear var x y ts_list

% disp(['Ad-hoc: Now duplicate the data'])
% TrajData.x = [TrajData.x, TrajData.x];
% TrajData.y = [TrajData.y, TrajData.y];

disp('Downsampled Trajectories: ')
% Note that x,y,ts_list clsoe are raw
[Nts, nparticles] = size(TrajData.x);
NJumps = (Nts-1)*nparticles;


disp(['NJumps = ', num2str(NJumps), '; Nts = ', num2str(Nts), '; nparticles = ', num2str(nparticles)]);

DownSampling_Param = struct('veloc_Profile', veloc_Profile, ...
                            'SamplingInterval_vis', SamplingInterval_vis, ...
                            'SamplingInterval', SamplingInterval, ...  % Obtained when Downsampling
                            't_offset', t_offset, ...
                            'Nts_RAW', Nts_RAW, 'nparticles_RAW', nparticles_RAW, ...
                            'NJumps_RAW', NJumps_RAW, ...
                            'Nts', Nts, 'nparticles', nparticles, ...
                            'NJumps', NJumps, 'DS_rate', DS_rate, ...
                            'traj_fullpath', traj_fullpath);
                        

%% Pre-process data
% Setup Rectangular Mesh
run('./Scripts/Script_SetupRectMesh');
disp('Finished Setting Up RectMesh.')
                        

%% Set up strings for storage of Pre-processed data
DataTypeString = 'BinnedTrajData';
run('./Scripts/Script_Filenames');

%% Position the trajectories w.r.t. Rectangular Grid
adhoc_off = 0;
adhoc_False = (1==0);
if strcmp(filename_BTD, filename_BTD_old) && adhoc_False
    %     load(filename_BTD)
    %     disp('Finished loading TrajJumps')
    disp('BinnedTrajData is already in memory. No need to sort it again.')
else
    disp('BinnedTrajData is NOT in memory. Now need to sort it.')
    if ~(adhoc_off == 1)
        [GridInd, alphax, alphay] = positioning_RectMesh(TrajData.x, TrajData.y, RectMesh_Param);
        
        TrajPosition = struct('GridInd', GridInd, ...
            'alphax', alphax, 'alphay', alphay );
        disp('Finished positioning and computing alpha.')
    else
        disp('TURN-OFF: positioning and computing alpha.')
    end
    
    %% Compute the transition matrix
    DataTypeString = 'BinnedTrajTransMat';
    run('./Scripts/Script_Filenames');

    
    if exist(filename_BTTM, 'file') == 2  && adhoc_False
        disp('BinnedTrajTransMat file already exist.');
        
        BinnedTrajTransMat = load(filename_BTTM);
        TransMat = BinnedTrajTransMat.TransMat;
        
        [TransDen, pi_Stat, RatioRemain, RatioRemainNeigh] = obtain_appended_TransMat(TransMat, Mesh_Struct);
        
        % Save the TransMat
        [status, msg, msgID] = mkdir(filepath);
        save(filename_BTTM, 'TransMat', 'TransDen', 'pi_Stat', ...
            'RatioRemain', 'RatioRemainNeigh', 'DownSampling_Param', ...
            'Mesh_Struct', 'traj_fullpath', '-v7.3')
        disp(['Appended  BinnedTrajTransMat - Saved a new file.'])

        
    else
        tic;
        [TransMat] = obtain_TransMat(TrajData, TrajPosition, Mesh);
        Dt_TransMat = toc;
        disp(['Finished Evaluating the Transition Matrix, under time = ', num2str(Dt_TransMat), 'sec'])
        
        [TransDen, pi_Stat, RatioRemain, RatioRemainNeigh] = obtain_appended_TransMat(TransMat, Mesh_Struct);
        disp(['Finished Appending the Transition Matrix.'])
        
        % Save the TransMat
        [status, msg, msgID] = mkdir(filepath);
        save(filename_BTTM, 'TransMat', 'TransDen', 'pi_Stat', ...
            'RatioRemain', 'RatioRemainNeigh', 'DownSampling_Param', ...
            'Mesh_Struct', 'traj_fullpath', '-v7.3')
        clear var TrajJumps_iBased
    end
    
    % Restore output to BinnedTrajData
    DataTypeString = 'BinnedTrajData';
    run('./Scripts/Script_Filenames');
    %% Label the jump to individual cells   
    % For Local approach only (assumed Ito interpretation: associate jumps to initial point)
    % Locate in which cell each particle lies in
    tic;
    %TrajJumps = sort_TrajJumps_Ito(TrajData, TrajPosition, Mesh);
    % Assumed uniform sampling interval: to reduce computation cost for Linear SDE
    disp('Assumed Uniform Time Step from Raw Data.')
    
    if ~(adhoc_off == 1)
        clear var TrajJumps NJumps_cell
        [TrajJumps, NJumps_cell] = sort_TrajJumps_Ito_MultipleSteps(TrajData, TrajPosition, Mesh, Nlevel);
        Dt_sortJumps= toc;
        assert(all(NJumps_cell >= 3))
        
        % For non-local approach: 21/9/2017
        % [TrajJumps, TrajJumps_Length] = sort_TrajJumps_iBased(TrajData, TrajPosition, Mesh);
        
        disp(['Finished sorting TrajJumps by Ito interpretation, under time = ', num2str(Dt_sortJumps), 'sec'])
        %%% Reminder: TrajJumps(cell_k) contains alphax, alphay, diffx, diffy, h
        
        if Nlevel == 1
            [TrajJumps_DA] = obtain_TrajJumps_Ito_DataAnalysis(TrajJumps, Mesh);
            disp('Finished Analyzing TrajJumps from Ito Statistically.')
            
            [status, msg, msgID] = mkdir(filepath);
            save(filename_BTD, 'TrajJumps_DA', 'DownSampling_Param', ...
                'Mesh_Struct', 'traj_fullpath', '-v7.3')
            disp('Finished saving TrajJumps.')
            
            clear var TrajJumps_DA
        end
    else
        disp(['Turned off sorting jumps.'])
    end
    
    filename_BTD_old = filename_BTD;    
end