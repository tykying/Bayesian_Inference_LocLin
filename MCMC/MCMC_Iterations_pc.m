%% MCMC Algorithm
%% Computation: G.1 MCMC move: burn-in
function [theta_Stat, MCMC_Stat] = MCMC_Iterations_pc(TrajJumps, theta_Stat, Prior_Param, Mesh, MHRW_Param, MCMC_Param)
Ncells = length(Mesh);

NVars = MCMC_Param.NVars;

theta_SSD = MHRW_Param.theta_SSD;

% Initialise theta
theta = theta_Stat.theta_Init;
theta_store = theta_Stat.theta_store;

% level of sampling
Nlevel = size(TrajJumps, 2);

Nsteps_pc = MCMC_Param.Nsteps_pcv * NVars;

Ntheta_store = theta_Stat.Ntheta_store;     % Include initial sample
assert(Ntheta_store == (Nsteps_pc+1));

% Reset the size of logPost_ts
% Compute the Initial log-Posterior, based on the given discretised fields
logPost_ts = zeros(Ncells, Ntheta_store);

acceptance_ratio = zeros(Ncells, NVars);

disp('Running MCMC_Iterations_pc');

% Clear memory
clear var TrajJumps_cellk TrajJumps_cellk_level

parfor cell_k = 1:Ncells
%for cell_k = 1:Ncells    
    TrajJumps_cellk = TrajJumps(cell_k, :);
    
    logPost_ts_cellk = logPost_ts(cell_k, :);

       
    % Need to initialise with fixed parameters for the ease of
    % parallelisation code to recognise it
    theta_SSD_cellk = zeros(1, NVars);
    theta_cellk = zeros(1, NVars);
    theta_prop_cellk = zeros(1, NVars);
    
    theta_SSD_cellk(1, :) = theta_SSD(cell_k, :);
    theta_cellk(1, :) = theta(cell_k, :);
    theta_prop_cellk(1, :) = theta(cell_k, :);

    theta_store_cellk = zeros(1, NVars, Ntheta_store);  % size = [1, 9, Ntheta_store]
    
    theta_store_cellk(1,:,:) = theta_store(cell_k, :, :);
    
    % Construct Proposal Probabilities
    logPrior_cellk = 0;
    logLik_cellk = 0;
    for level = 1:Nlevel
        TrajJumps_cellk_level = TrajJumps_cellk(level);
               
        logLik_cellk = logLik_cellk + compute_logLik_jump_cell_linearSDE(TrajJumps_cellk_level, theta_cellk);
    end
    
    logPrior_cellk = compute_logPrior_cell(Prior_Param, theta_cellk);
    logPost_cellk = logPrior_cellk + logLik_cellk;
    
    logPost_ts_cellk(1) = logPost_cellk;
    
    % Begin the MCMC-MH
    acceptance_counter = zeros(1, NVars);
    rejection_counter = zeros(1, NVars);

    
    for step_counter = 1:Nsteps_pc       
        % MHRW: Proposal
        % All parameters at a time
        %theta_prop_cell = theta(cell_k, :) + theta_SSD(cell_k, :).*randn(1, NVars);  % propose new U, K
        
        % One parameter at a time
        cmpn = mod(step_counter-1, NVars)+1;
        theta_prop_cellk(1, cmpn) = theta_cellk(1, cmpn) + theta_SSD_cellk(1, cmpn).*randn(1, 1);  % propose new U

        % Ensure Posivities
        theta_prop_cellk = restrict_theta_cell(theta_prop_cellk);
        
        % Compute the log-Post of each cell
        % Idea: Cell-based summation (For non-U1K0: consider support of basis)
        % Localised -> Prior and Likelihood does not depends on other cells
        logLik_prop_cellk = 0;
        for level = 1:Nlevel
            TrajJumps_cellk_level = TrajJumps_cellk(level);

            logLik_prop_cellk = logLik_prop_cellk + compute_logLik_jump_cell_linearSDE(TrajJumps_cellk_level, theta_prop_cellk);
        end
        logPrior_prop_cellk = compute_logPrior_cell(Prior_Param, theta_prop_cellk);
        logPost_prop_cellk = logPrior_prop_cellk + logLik_prop_cellk;
        
        % Accept a new theta
        uniformRV = rand;    % uniform random variable (c.f.  p.115 J.S Liu)
        
        q21 = 1; q12 = 1;
        alphaUniRV = min(1, exp(logPost_prop_cellk-logPost_cellk)*q21/q12);   % P(uniformRV < alphaUniRV) = alphaUniRV; accept-reject step
        
        if (uniformRV < alphaUniRV)
            % Update the grid values of U, V, K
            theta_cellk = theta_prop_cellk;
            
            logPost_cellk = logPost_prop_cellk;
            
            acceptance_counter(cmpn) = acceptance_counter(cmpn) + 1;
        else
            rejection_counter(cmpn) = rejection_counter(cmpn) + 1;
        end
        
       % Storage of theta and logPost
        logPost_ts_cellk(step_counter+1) = logPost_cellk;
        theta_store_cellk(1,:,step_counter+1) = theta_cellk;        
    end
    
    % Export back to global variable
    logPost_ts(cell_k, :) = logPost_ts_cellk;    
    theta_store(cell_k, :, :) = theta_store_cellk;  % size = [1, 9, N_eff_sample]
    
    % Compute acceptance ratio; Export to global variable
    %disp(['AccCounter, RejCounter = ', num2str([acceptance_counter, rejection_counter]));
    acceptance_ratio(cell_k, :) = acceptance_counter./(acceptance_counter+rejection_counter);
    
    
    if (mod(cell_k, 64) == 1)
        disp(['MCMC: ', num2str(Nsteps_pc), ' steps done in cell ', num2str(cell_k)]);
    end
end

% Output
MCMC_Stat = struct('acceptance_ratio', acceptance_ratio, ...
    'Ntheta_store', Ntheta_store, ...
    'version', 'pc_parallel');

theta_Stat.theta_store = theta_store;
theta_Stat.logPost_ts = logPost_ts;
end
