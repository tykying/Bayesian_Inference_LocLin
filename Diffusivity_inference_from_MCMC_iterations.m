%% Inference
% Bayesian estimates
% Initialise the MCMC Algorithm
NVars = 9;   % Hardcoded
Nsteps = Nsteps_pcv*length(Mesh)*NVars;

% Specify prior
% Reminder: variance = sqrt(kappa)^2
Prior_Param = struct(...
    'b', [0, 1], ...
    'nu', [-1E-5, 1E-5], ...
    'tr', [-1E-5, 1E-5], ...
    'sigma', [sqrt(10), sqrt(40000)], ...
    'PriorType', 2);  % 2: Uniform prior: from [sqrt(10), 200] for sigmas

MCMC_Param = struct('InfScheme', InfScheme, ...
    'Nsteps_pcv', Nsteps_pcv, ...
    'NJumps_cell', NJumps_cell, ...
    'space_scale', space_scale, ...
    'time_scale', time_scale, ...
    'Pe', Pe, ...
    'NVars', NVars, ...
    'veloc_Profile', veloc_Profile);

%% Determine MHRW_Param and theta_Init
adhoc_False = (1==0);

DataTypeString = 'BayesMHConfig';
run('./Scripts/Script_Filenames');
if (exist(filename_BMConfig, 'file') == 2) && adhoc_False
    % Optimal theta_SSD already determined; Update MHRW_Param
    load(filename_BMConfig);
else
    disp(['Spinning Up MCMH']); 
    
    NTrials = length(nTrial_List);
    theta_SSD_List = zeros(length(Mesh), NVars, NTrials);
    theta_Init_List = zeros(length(Mesh), NVars, NTrials);
    AccRatio_List = zeros(length(Mesh), NVars, NTrials);
    logPost_MAP_List = zeros(length(Mesh), NTrials);
        
    
    % Spin up with 2^16 steps
    MCMC_Param_SpinUp = MCMC_Param;
    %MCMC_Param_SpinUp.Nsteps_pc = 2^18;
    
    % Fine Tuning
    MCMC_Param_TUNE = MCMC_Param;
    MCMC_Param_TUNE.Nsteps_pcv = 1024;
    
    for nTrial = nTrial_List
        rng(nTrial*48);
        
        % Deterministic initialisation
        [theta_Init, theta_Stat, MHRW_Param] = initialise_MCMC(MCMC_Param_SpinUp, Mesh_Struct, TrajJumps);
        
        DataTypeString = 'BayesMHRW';
        run('./Scripts/Script_Filenames');
        
        if (exist(filename_BMHRW) == 2) && adhoc_False
            load(filename_BMHRW);
        else
            tic;
            [theta_Stat, MCMC_Stat] = MCMC_Iterations_pc(TrajJumps, theta_Stat, Prior_Param, Mesh, MHRW_Param, MCMC_Param_SpinUp );
            Dt_Sampling = toc;
            disp(['Finished MCMC Spin Up, under time = ', num2str(Dt_Sampling), 'sec']);
            
            % MAP Per sessions
            % First entry must be 1 for the Global;
            % Second entry can be any positive
            NSessions = 64;
            [theta_MAP_Sessions, logPost_MAP_Sessions] = locate_theta_MAP_Sessions(theta_Stat.theta_store, theta_Stat.logPost_ts, NSessions);
            [theta_MAP_Global, logPost_MAP_Global] = locate_theta_MAP_Sessions(theta_Stat.theta_store, theta_Stat.logPost_ts, 1);
            
            % Tune for theta_MAP
            theta_Init_TUNE = theta_MAP_Global;
           
            [theta_Init_TUNE, theta_Stat, MHRW_Param] = initialise_MCMC(MCMC_Param_TUNE, Mesh_Struct, TrajJumps, theta_Init_TUNE, MHRW_Param.theta_SSD);
            
            tic;
            [theta_Stat_TUNE, MHRW_Param_TUNE, MCMC_Stat_TUNE] = tune_MHRW_Param(TrajJumps, theta_Init_TUNE, Prior_Param, Mesh_Struct, MHRW_Param, MCMC_Param_TUNE, MCMC_Stat);
            Dt_Sampling = toc;
            disp(['Finished tuning MHRW_Param, under time = ', num2str(Dt_Sampling), 'sec']);
            
            % [min(AccRatio_I{end},[], 2), max(AccRatio_I{end},[], 2)]
            
            [status,msg,msgID] = mkdir(filepath);
            save(filename_BMHRW, 'theta_Stat_TUNE', 'MHRW_Param_TUNE', 'MCMC_Stat_TUNE', ...
                'theta_MAP_Global', 'logPost_MAP_Global', ...
                'theta_MAP_Sessions', 'logPost_MAP_Sessions', ...
                'theta_Init', 'MHRW_Param', ...
                'MCMC_Param_SpinUp', 'Prior_Param', 'DownSampling_Param', ...
                'Mesh_Struct', 'traj_fullpath', '-v7.3');
        end
                
        fprintf('Trial %d: Obtained initialised theta_Stat. \n', nTrial)
        
        % Record the final iteration
        if ~isempty(MHRW_Param_TUNE)
            theta_SSD_MAP = MHRW_Param_TUNE{end}.theta_SSD;
            AccRatio_MAP = MCMC_Stat_TUNE{end}.acceptance_ratio;
        else
            theta_SSD_MAP = MHRW_Param.theta_SSD;
            AccRatio_MAP = MCMC_Stat.acceptance_ratio;
        end
        
        theta_SSD_List(:, :, nTrial) = theta_SSD_MAP;
        theta_Init_List(:, :, nTrial) = theta_MAP_Global;
        AccRatio_List(:, :, nTrial) = AccRatio_MAP;
        logPost_MAP_List(:, nTrial) = logPost_MAP_Global;
    end
        
    %% Pick the one with maximal logPost 
    % M1: Use MAP estimates
    %[logPost_MAP_Picked, nTrial_Picked] = max(logPost_MAP_List, [], 2);
    [logPost_MAP_Sorted, nTrial_MAP_Sorted] = sort(logPost_MAP_List, 2, 'descend');

    % M2: Use smallest possible U, V
    % Compute b magnitude
    b_Mag_List = zeros(length(Mesh), length(nTrial_List));
    for nTrial = 1:length(nTrial_List)
        b_Mag = theta_Init_List(:, 1, nTrial);
        b_Mag_List(:, nTrial) = squeeze(b_Mag);
    end
    
    [b_Mag_List_Sorted, nTrial_Sorted] = sort(b_Mag_List, 2, 'ascend');
    
    theta_Init_Picked = zeros(length(Mesh), NVars);
    theta_SSD_Picked = zeros(length(Mesh), NVars);
    AccRatio_Picked = zeros(length(Mesh), NVars);
    nTrial_Picked = zeros(length(Mesh), 1);
    for cell_k = 1:length(Mesh)
%         nTrial_ind = 1;
         
         nTrial_Picked_cell_k = nTrial_Sorted(cell_k, 1);
%         AccRatio_cellk_UV_Picked = squeeze(AccRatio_List(cell_k, 1:2, nTrial_Picked_cell_k));
% 
%         while any(AccRatio_cellk_UV_Picked > 0.4) || any(AccRatio_cellk_UV_Picked < 0.1)
%             nTrial_ind = nTrial_ind + 1;
%             
%             nTrial_Picked_cell_k = nTrial_Sorted(cell_k, nTrial_ind);
%             AccRatio_cellk_UV_Picked = squeeze(AccRatio_List(cell_k, 1:2, nTrial_Picked_cell_k));
%         end
        
        nTrial_MAP_cell_k = nTrial_MAP_Sorted(cell_k, 1);

         if (nTrial_Picked_cell_k ~= nTrial_MAP_cell_k)
             fprintf('cell %d: nTrial MAP = %d <> LowUV = %d. \n', cell_k, nTrial_MAP_cell_k, nTrial_Picked_cell_k);
             % Picked nTrial_MAP anyway
             nTrial_Picked_cell_k = nTrial_MAP_cell_k;
             fprintf('[Picked nTrial_MAP] b_Mag_List =  ');
             disp(b_Mag_List(cell_k, :));
         end
        
        nTrial_Picked(cell_k) = nTrial_Picked_cell_k;
        theta_Init_Picked(cell_k, :) = squeeze(theta_Init_List(cell_k, :, nTrial_Picked_cell_k));
        theta_SSD_Picked(cell_k, :) = squeeze(theta_SSD_List(cell_k, :, nTrial_Picked_cell_k));
        AccRatio_Picked(cell_k, :) = squeeze(AccRatio_List(cell_k, :, nTrial_Picked_cell_k));
    end

    theta_Init = theta_Init_Picked;
    MHRW_Param.theta_SSD = theta_SSD_Picked;    

    DataTypeString = 'BayesMHConfig';
    run('./Scripts/Script_Filenames');
    [status,msg,msgID] = mkdir(filepath);
    save(filename_BMConfig, 'MHRW_Param', 'theta_Init', ...
        'nTrial_Picked', 'theta_Init_Picked', 'theta_SSD_Picked', ...
        'AccRatio_Picked', 'logPost_MAP_List', ...
        'theta_SSD_List', 'theta_Init_List', 'AccRatio_List', ...
        'MCMC_Param_SpinUp', 'MCMC_Param_TUNE', ...
        'Prior_Param', ...
        'DownSampling_Param', ...
        'Mesh_Struct', 'traj_fullpath', '-v7.3')
end

%% Bayesian Inference
DataTypeString = 'BayesSampleStat';
run('./Scripts/Script_Filenames');

for nTrial = nTrial_List
    %% Iterations
    if (exist(filename_BSS, 'file') == 2) && adhoc_False
        disp(['BayesSampleStat already exist. Now Loading.']);
        load(filename_BSS);
    else
        rng(nTrial*611);
        
        % Choose INIT
        theta_Init = theta_Init_List(:, :, nTrial);
        theta_SSD = theta_SSD_List(:, :, nTrial);
        
        [theta_Init, theta_Stat, MHRW_Param] = initialise_MCMC(MCMC_Param, Mesh_Struct, TrajJumps, theta_Init, theta_SSD);
        
        tic;
        [theta_Stat, MCMC_Stat] = MCMC_Iterations_pc(TrajJumps, theta_Stat, Prior_Param, Mesh, MHRW_Param, MCMC_Param );
        Dt_Sampling = toc;
        disp(['Finished MCMC Sampling of trial ', num2str(nTrial), ', under time = ', num2str(Dt_Sampling), 'sec']);
        
        % Mark the theta_MAP
        [theta_MAP, logPost_MAP, MAP_SampleInd] = locate_theta_MAP(theta_Stat.theta_store, theta_Stat.logPost_ts);
        theta_Stat.theta_MAP = theta_MAP;
        theta_Stat.logPost_MAP = logPost_MAP;
        
        % Compress theta_store to save storage space
        theta_Stat.theta_store  = zip_theta(theta_Stat.theta_store);
        
        [status,msg,msgID] = mkdir(filepath);
        save(filename_BSS, 'theta_Stat', 'MCMC_Stat', 'MHRW_Param', ...
            'MCMC_Param', 'theta_Init', 'Prior_Param', 'DownSampling_Param', ...
            'Mesh_Struct', 'traj_fullpath', '-v7.3')
    end
    
    %% Data Storage in Memory
    if (store_expt == 1)
        theta_Stat_List{end+1} = theta_Stat;
        DownSampling_Param_List{end+1} = DownSampling_Param;
        MCMC_Stat_List{end+1} = MCMC_Stat;
        MCMC_Param_List{end+1} = MCMC_Param;
    else
        clear var theta_Stat
    end
end