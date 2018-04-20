%% Initialise the MCMC algorithm
% Specifications:

% Input:
% Nx_cell, Ny_cell: Number of cells in Spatial Discretisation
% kappa_scale, u_scale, v_scale: For initialisation and setup prior
% length(Jumps): For initialisation proposal variables

% Output:
% Jcell_t0(i): A vector that tells the cell (k) where the jump(i) comes from

%% Computation
function [theta_Init, theta_Stat, MHRW_Param] = initialise_MCMC(MCMC_Param, Mesh_Struct, TrajJumps, theta_Init_IN, theta_SSD_IN)
% Initialisation
Mesh = Mesh_Struct.Mesh;
Nsteps_pcv = MCMC_Param.Nsteps_pcv;
InfScheme = MCMC_Param.InfScheme;
NVars = MCMC_Param.NVars;
space_scale = MCMC_Param.space_scale;
time_scale = MCMC_Param.time_scale;
veloc_Profile = MCMC_Param.veloc_Profile;


Ncells = length(Mesh);
Nsteps_pc = Nsteps_pcv*NVars;

% space, time: SI unit
theta_Scale = zeros(Ncells, NVars);
theta_Scale(:, 1) = 0.01 *(time_scale/space_scale);
theta_Scale(:, [3:4,6]) = 1E-7 *(time_scale);
theta_Scale(:, 7:8) = 10 *sqrt(time_scale/space_scale^2);
theta_Scale(:, [2,5,9]) = 1;  % Angles

if nargin == 3    
    zero_init = 1;
    if zero_init == 1
        disp('theta_init: Zero Initial Condition.');
       
        % Set to be zeros
        theta_Init = zeros(Ncells, NVars);

        theta_Init(:, 7) = 50 *sqrt(time_scale/space_scale^2);
        theta_Init(:, 8) = 30 *sqrt(time_scale/space_scale^2);
        
        theta_SSD = theta_Scale*0.1;
    else
        disp('initialise_MCMC: theta_init given by Linear Fit.');
        
        % theta_init: struct; contains all discrete fields globally
        [theta_Init, theta_Scale] = initialise_theta_pwl_randomised(MCMC_Param, Mesh, TrajJumps);
        
        % Compute the standard deviation for the RW samplers
        % N = 4*NCells; I = (2*u_scale)^(-3/2); 
        % l = 2.38/sqrt(2*I); U_sampling_sd = sqrt(2/N)*l;
        
        % Allow each cell to have its own sampling sd
        theta_SSD = abs(theta_Scale)*0.05;
        theta_SSD(:,9) = pi/72;
    end
    
elseif nargin == 5
    disp('theta_init and theta_SSD given by Input.');
    
    theta_Init = theta_Init_IN;
    theta_SSD = theta_SSD_IN;
end


% InfScheme
if contains(InfScheme, 'PWC')
    theta_SSD(:, 3:6) = 0;
    theta_Init(:, [3:6]) = 0;
elseif contains(InfScheme, 'INC')
    theta_SSD(:, 6) = 0;
    theta_Init(:, 6) = 0;
end

%% Output:
% Storage of samples
Ntheta_store = Nsteps_pc+1;  % +1: Store initial condition

% Set-up theta_store: Initialise it in iterations
theta_store = zeros([size(theta_Init), Ntheta_store]);
theta_store(:,:,1) = theta_Init;

theta_Stat = struct( ...
    'theta_Init', theta_Init, ...
    'Nsteps_pcv', Nsteps_pcv, ...
    'Nsteps_pc', Nsteps_pc, ...
    'Ntheta_store', Ntheta_store, ...
    'theta_Scale', theta_Scale, ...
    'theta_store', theta_store);

MHRW_Param = struct('theta_SSD', theta_SSD);

end