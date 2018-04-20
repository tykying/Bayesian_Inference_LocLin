%clear all
close all
%% Initialisation
traj_fullpath_old = '';
filename_BTD_old = '';

%% Parameters Available
% Non-QGM2
veloc_Profile_List = {'shear', 'linear_shear', 'const_skew', 'childress_soward', 'taylor_green_noisy', 'childress_soward'};

% QGM2
veloc_Profile_List = {'QGM2_DSpart', 'QGM2_DSmSep',  ...        % Long Trajectories
                    'QGM2_nReal_DSSpat', 'QGM2_nReal_DSTemp'};  % Multiple Realisations

                
InfScheme_List = {'LocLin', 'LocLinINC', 'LocLinPWC'};
RunProfile_List = {'linear_jet', 'QGM2_DStemp', 'linear_shear', 'cosine_TC_pwc', 'cosine_TC', 'tg_w_mean_TC_pwc_vs_SIv', 'tg_w_mean_TC_pwc_vs_Npart'};
%% Setup the RunProfile
RunProfile ='QGM2_DStemp';
%RunProfile ='linear_jet';

run('./Scripts/Script_Parameters.m');

% SamplingInterval_vis_List = [64];
% nTrial_List = [1:3];
% store_expt = 1;
% InfScheme_List = {'LocLinINC'};
% veloc_Profile_List = {'linear_jet_CplxEigen_InC_ISO'};
% Nsteps_pcv_List = [2048*4*8];

%% Data storage
theta_Stat_List = {};
DownSampling_Param_List = {};
MCMC_Stat_List = {};
MCMC_Param_List = {};

%% Main Loop for Diffusivity and Binned Trajectories Transition Matrix
for u_scale_i = 1:length(u_scale_List)  % ad-hoc for TG with drift
    u_scale = u_scale_List(u_scale_i);
    
for u_angle = u_angle_List
    disp(['[U, phi] = ', num2str([u_scale, u_angle])])
    
for Layer = Layer_List
for veloc_Profile_iterator = 1:numel(veloc_Profile_List)
% From below onwards: based on same trajectories
for SamplingInterval_vis = SamplingInterval_vis_List
for DS_rate = DS_rate_List
for Nx_cell_ARG = Nx_cell_ARG_List
for InfScheme_iterator = 1:numel(InfScheme_List)  % Exact Inference vs Euler Approximation; TODO: Can be move to innest loop
for Nsteps_pcv = Nsteps_pcv_List  % Number of MCMC Samples per cell

veloc_Profile = veloc_Profile_List{veloc_Profile_iterator};
InfScheme = InfScheme_List{InfScheme_iterator};

run('Diffusivity_inference_from_MCMC_PreProcess');
run('Diffusivity_inference_from_MCMC_iterations');


end
end
end
end
end
end
end
end
end

%% Make Use of Stored Data
if store_expt == 1    
    %run('./Scripts/Script_OutputFigures.m');
end

for nTrial = 1:3
    run('./data_visualisation/prelim_check_BayesSampleStat_storedSamples.m');
end
