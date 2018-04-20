% % % Good example for demonstration: 
% % % Domain considered: L = 4;
% % % 1) childress_soward, sinusoidal; Pe=10; SIv=32; Nx_cell=16; => show how crude method fails
% % % 2) taylor_green_noisy, const; Pe=10; SIv=32; Nx_cell=16; => Able to capture gradient when sufficient data?
% % 
% % % Mix Trajectories:
% % % 1) + 2) Mix together => Does not seem to make the diffusivity diverge; Probably because after all the jumps are of the same scale

%% Parameters Available
% Non-QGM2
veloc_Profile_List = {'shear', 'linear_shear', 'const_skew', 'childress_soward', 'taylor_green_noisy', 'childress_soward'};

% QGM2
veloc_Profile_List = {'QGM2_DSpart', 'QGM2_DSmSep',  ...        % Long Trajectories
                    'QGM2_nReal_DSSpat', 'QGM2_nReal_DSTemp'};  % Multiple Realisations

                
InfScheme_List = {'LocLin', 'LocLinINC', 'LocLinPWC'};

%% Begin
% Number of Initial Conditions
nTrial_List = [1];

% Default value
u_angle_List = [0];
u_scale_List = [1];


if (strcmp(RunProfile, 'QGM2_DStemp'))
    veloc_Profile_List = {'QGM2_DStemp'};
    kappa_Profile = 'const';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters
    % MCMC Parameter
    InfScheme_List = {'LocLinINC'};
    Nsteps_pcv_List = [1024*(2^5)];  % Included the burnt-in samples
    Nsteps_pcv_List = [50000];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [16];
    DS_rate_List = [1];

    
    % Key Parameter
    %SamplingInterval_vis_List = [1, 16, 32:32:256];
    nTrial_List = [1:5];
    
    % Store theta output
    store_samples = 1;
    store_expt = 0;
end


if (strcmp(RunProfile, 'linear_shear'))
    veloc_Profile_List = {'linear_shear'};
    kappa_Profile = 'const_skew';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters
    % MCMC Parameter
    InfScheme_List = {'Local'};
    %Nsteps_pcv_List = [1024*(2^9)];  % Included the burnt-in samples
    Nsteps_pcv_List = [1024*(2^6)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [1];
    DS_rate_List = [64];

    
    % Key Parameter
    %SamplingInterval_vis_List = [256:-32:32, 16];
    %SamplingInterval_vis_List = [96];
    SamplingInterval_vis_List = [16];
    nTrial_List = [1:5];
    
    % Store theta output
    store_samples = 1;
    store_expt = 1;
end

if (contains(RunProfile, 'linear_jet'))
    veloc_Profile_List = {'linear_jet_RealEigen', 'linear_jet_CplxEigen'};
    veloc_Profile_List = {'linear_jet_RealEigen_InC', 'linear_jet_CplxEigen_InC'};
    veloc_Profile_List = {'linear_jet_RealEigen_ISO', 'linear_jet_CplxEigen_ISO'};

    kappa_Profile = 'const_skew';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters
    % MCMC Parameter
    InfScheme_List = {'LocLin'};
    Nsteps_pcv_List = [1024*(2^6)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [1];
    DS_rate_List = [1];

    
    % Key Parameter
    %SamplingInterval_vis_List = [256:-32:32, 16];
    nTrial_List = [1:1];
    
    % Store theta output
    store_samples = 1;
    store_expt = 1;
end

if (strcmp(RunProfile, 'cosine_TC')  || strcmp(RunProfile, 'cosine_TC_pwc'))
    veloc_Profile_List = {'cosine'};
    kappa_Profile = 'const';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters
    % MCMC Parameter
    InfScheme_List = {'Local'};
    Nsteps_pcv_List = [1024*(2^7)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [16];
    DS_rate_List = [1];
    nTrial_List = [1];

    
    % Key Parameter
    SamplingInterval_vis_List = [1,2, 4:4:512, 520:8:1024];
    SamplingInterval_vis_List = 2.^[0:10]; %[1,2, 4:4:512, 520:8:1024];

    %[32:32:1024, 1024:512:4096];  %% To illustrate importance of resolving shear
    
    % Store theta output
    store_samples = 1;
    store_expt = 1;
end


if strcmp(RunProfile, 'childress_soward_TC')   
    veloc_Profile_List = {'childress_soward'};
    kappa_Profile = 'sinusoidal';    % Hardcoded
    Pe = 10;     % Hardcoded
    
    % Parameters
    % MCMC Parameter
    InfScheme_List = {'Local'};
    Nsteps_pcv_List = [1024*(2^4)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [16];
    DS_rate_List = [1];

    
    % Key Parameter
    SamplingInterval_vis_List = [64];  %% To illustrate importance of resolving shear
    
    % Store theta output
    store_samples = 1;
    store_expt = 1;
end

if strcmp(RunProfile, 'childress_soward_TC_vs_shear')   
    veloc_Profile_List = {'childress_soward'};
    kappa_Profile = 'sinusoidal';    % Hardcoded
    Pe = 10;     % Hardcoded
    
    % Parameters
    % MCMC Parameter
    InfScheme_List = {'Local'};
    Nsteps_pcv_List = [1024*(2^4)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [16];
    DS_rate_List = [1];


    % Key Parameter
    SamplingInterval_vis_List = [32];  %% To illustrate importance of resolving shear
    
    
    % Store theta output
    store_samples = 0;
    store_expt = 1;
end


if strcmp(RunProfile, 'tg_w_mean_TC_TEST')   
    veloc_Profile_List = {'tg_w_mean'};
    kappa_Profile = 'const';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters 
    u_angle_List = linspace(0, pi/2, 16+1); %!KEY!%
    u_angle_List = u_angle_List(7);   %!KEY!%  % Only choose the value 6pi/32
    u_angle_List = 0;
    
    u_scale_List = [0];   %!KEY!%

    % MCMC Parameter
    InfScheme_List = {'Local'};
    Nsteps_pcv_List = [1024*(2^7)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [1];
    DS_rate_List = [1];   %!KEY!%
    SamplingInterval_vis_List = [1024];  %% To illustrate importance of resolving shear
    
    
    % Store theta output
    store_samples = 1;
    store_expt = 1;
end



if (strcmp(RunProfile, 'tg_w_mean_TC_vs_SIv') || strcmp(RunProfile, 'tg_w_mean_TC_pwc_vs_SIv'))
    veloc_Profile_List = {'tg_w_mean'};
    kappa_Profile = 'const';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters 
    u_angle_List = linspace(0, pi/2, 16+1); %!KEY!%
    u_angle_List = u_angle_List(7);   %!KEY!%  % Only choose the value 6pi/32

    u_scale_List = [1];   %!KEY!%

    % MCMC Parameter
    InfScheme_List = {'Local'};

    Nsteps_pcv_List = [1024*(2^7)];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [1];
    DS_rate_List = [4];   %!KEY!%
    SamplingInterval_vis_List = [1:12, 14:2:32, 64, 96, 128, 192, 256, 384, 512, 768, 1024];  %% To illustrate importance of resolving shear
    SamplingInterval_vis_List = [8:8:2048];  %% To illustrate importance of resolving shear

    % Store theta output
    store_samples = 1;
    store_expt = 1;
end

if (strcmp(RunProfile, 'tg_w_mean_TC_vs_Npart') || strcmp(RunProfile, 'tg_w_mean_TC_pwc_vs_Npart'))
    veloc_Profile_List = {'tg_w_mean'};
    kappa_Profile = 'const';    % Hardcoded
    Pe = 100;     % Hardcoded
    
    % Parameters 
    u_angle_List = linspace(0, pi/2, 16+1); %!KEY!%
    u_angle_List = u_angle_List(7);   %!KEY!%  % Only choose the value 6pi/32

    u_scale_List = [1];   %!KEY!%

    % MCMC Parameter
    InfScheme_List = {'Local'};
    Nsteps_pcv_List = [1024*1024];  % Included the burnt-in samples

    % Downsampling
    Layer_List = [2];   % Irrelevant

    Nlevel = 1;  % Irrelevant
    Nx_cell_ARG_List = [1];
    DS_rate_List = [64, 16, 4];   %!KEY!%
    SamplingInterval_vis_List = [32/(2^-5)];  %h=2^-5
    
    % Store theta output
    store_samples = 1;
    store_expt = 1;
end
