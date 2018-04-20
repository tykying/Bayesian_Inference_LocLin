function [theta_Stat_TUNE, MHRW_Param_TUNE, MCMC_Stat_TUNE] = tune_MHRW_Param(TrajJumps, theta_Init, Prior_Param, Mesh_Struct, MHRW_Param, MCMC_Param, MCMC_Stat)
disp('Iterating for optimal MHRW theta_SSD')

Mesh = Mesh_Struct.Mesh;

MCMC_Param_const= MCMC_Param;

% Initialisation
MHRW_Param_i = MHRW_Param;
MCMC_Stat_i = MCMC_Stat;

theta_SSD_i = MHRW_Param.theta_SSD;
AccRatio_i = MCMC_Stat.acceptance_ratio;

% List of components considered
if contains(MCMC_Param_const.InfScheme, 'PWC')
    cmp_list = [1, 7:8];
elseif contains(MCMC_Param_const.InfScheme, 'INC')
    cmp_list = [1,3:4,6:8];
else
    cmp_list = [1:MCMC_Param_const.NVars];
end

    
[theta_Init_i, theta_Stat_i, MHRW_Param_i] = initialise_MCMC(MCMC_Param_const, Mesh_Struct, TrajJumps, theta_Init, theta_SSD_i);


AccRatio_MIN = 0.15;
AccRatio_MAX = 0.35;

AccRatio_i_F = AccRatio_i(:, cmp_list);  % F: Filtered
AccRatio_BLB = (AccRatio_i_F < AccRatio_MIN); % BLB: Below Lower Bound
AccRatio_AUB = (AccRatio_i_F > AccRatio_MAX); % AUB: Above Upper Bound

theta_Stat_TUNE = {};
MHRW_Param_TUNE = {};
MCMC_Stat_TUNE = {};

iteration = 0;
while ( any(AccRatio_BLB(:)) || any(AccRatio_AUB(:)) ) && (iteration < 16)
    rng(iteration*1991);

    iteration = iteration + 1;
    disp(['Tune MHRW iteration: ', num2str(iteration)])
        
    % Reiterate to update theta_SSD        
    theta_SSD = MHRW_Param_i.theta_SSD;
    theta_SSD(AccRatio_i < AccRatio_MIN) = theta_SSD(AccRatio_i < AccRatio_MIN) * (2/3);
    theta_SSD(AccRatio_i > AccRatio_MAX) = theta_SSD(AccRatio_i > AccRatio_MAX) * (4/3);
    theta_SSD_Sigmas = 0.5*(theta_SSD(:, 7) + theta_SSD(:, 8));
    theta_SSD(:, 7) = theta_SSD_Sigmas;
    theta_SSD(:, 8) = theta_SSD_Sigmas;
    
    
    MHRW_Param_i.theta_SSD = theta_SSD;
    
    [theta_Stat_i, MCMC_Stat_i] = MCMC_Iterations_pc(TrajJumps, theta_Stat_i, Prior_Param, Mesh, MHRW_Param_i, MCMC_Param_const );

    
    AccRatio_i = MCMC_Stat_i.acceptance_ratio;
    
    % Exclude angles in the decision
    AccRatio_i_F = AccRatio_i(:, cmp_list);
    AccRatio_BLB = (AccRatio_i_F < AccRatio_MIN); 
    AccRatio_AUB = (AccRatio_i_F > AccRatio_MAX);
    fprintf('Exclude Phi: Acceptance Ratio: [min, max] = [%f, %f]. \n', min(AccRatio_i_F(:)), max(AccRatio_i_F(:)) );
    
    theta_Stat_TUNE{end+1} = theta_Stat_i;
    MHRW_Param_TUNE{end+1} = MHRW_Param_i;
    MCMC_Stat_TUNE{end+1} = MCMC_Stat_i;
end

end