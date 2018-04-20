% Generate Simulated Trajectories
%% Estimation of a variable diffusivity from a Brownian trajectory

% RU1E OF THUMB: 
%  (t, dim, *sample*)
% 1) Increase in time = Increase in row 
% 2) If stationary, number of row <=> along x,y,z

% 3) Gradient of one scalar is ALWAYS row vectors at first. If it involves mU1tiple sampes then treated in case by case basis.

clear all
close all

rng(0)
Use_SavedSimTraj = 0;

%% Ideal testcase
RunProfile_List = {'linear_jet_RealEigen_ISO', 'linear_jet_CplxEigen_ISO', ...
                    'linear_jet_RealEigen_InC_ISO', 'linear_jet_CplxEigen_InC_ISO', ...
                    'linear_jet_RealEigen', 'linear_jet_CplxEigen', ...
                    'linear_jet_InC', 'cosine', 'tg_w_mean'};


RunProfile = 'linear_jet_CplxEigen_InC_ISO';
if contains(RunProfile, 'linear_jet')
    veloc_Profile = RunProfile;
    kappa_Profile = 'const_skew';
    
    L = 1E6;  % Bin size = L/10;
    GradU_scale = [1E-7, 1E-7];
    vortix_strength = 1;  % Irrelevant
    
    % Cell-centre velocity: 
    U_scale = 1E-1;
    b = U_scale*[3, 1];
    
    b_mag_List = [norm(b)];
    b_angle_List= [atan(b(2)/b(1))];
    
    Pe_List = [100];
    %Pe_List = [0.001, 0.1, 1, 10, 100, 1000];

    
    % Time stepping
    t_0 = 0; t_1 = 32*24*3600;
    npart_pc = 2^5;
    
    Nts= 2^10;
    h = (t_1-t_0)/Nts;
    sub_timestep = 16;
end


%% Initial conditions of particles
PART_ranges_min = [-L, -L]/10;
PART_ranges_max = [L, L]/10;

npart_pc_DIR = [npart_pc, npart_pc];

ncell_DIR = [1, 1];

%% Set-up variable names for output
Npart = prod(npart_pc_DIR.*ncell_DIR);
TempIntv = (t_1-t_0);

TempIntv_vis = (TempIntv/(24*3600));
SamplingInterval_vis = round(h/(24*3600));  % Could be rounded to 0
nReal = 1;   % irrelevant


assert(vortix_strength == 1);

% param_List = struct('kappa_scale', 0, 'L', L, 'Pe', Pe, 'U_scale', U_scale, ...
%     'U_angle', U_angle, 'GradU_scale', GradU_scale, 'vortix_strength', vortix_strength, ...
%     'Npart', Npart, 'TempIntv', TempIntv, ...
%     'SamplingInterval_vis', SamplingInterval_vis, 'nReal', nReal);

%% Experiment parameters
% Configure here to run for different parameters
k=0;
for Pe = Pe_List
    for b_mag = b_mag_List
        for b_angle = b_angle_List
            k = k+1;
            
            U_scale = 10^floor(log10(abs(b_mag)));
            
            kappa_scale = U_scale*L/Pe;
            sigma_scale = sqrt(2*kappa_scale);
            
            param_List(k) = struct('kappa_scale', kappa_scale, 'L', L, ...
                'Pe', Pe, 'b_mag', b_mag, 'b_angle', b_angle, ...
                'U_scale', U_scale, 'GradU_scale', GradU_scale, ...
                'vortix_strength', vortix_strength, ...
                'Npart', Npart, ...
                'TempIntv', TempIntv, 'TempIntv_vis', TempIntv_vis, ...
                'h', h, 'SamplingInterval_vis', SamplingInterval_vis, ...
                'nReal', nReal); 
        end
    end
end

%% Runs of experiments
trajdata_List(length(param_List)) = struct();

%parfor k = 1:length(param_List)
for k = 1:length(param_List)
param = param_List(k);

fprintf('U = %f \n', norm(param.U_scale, 2));

kappa_scale = param.kappa_scale;
sigma_scale = sqrt(2*kappa_scale);
% Reminder:
% kappa: diffusivity for FK
% sigma: diffusion coeff in SDE

kappa_fldStruct = Prescribed_DiffusivityTensorField(kappa_Profile, param);
veloc_fldStruct = Prescribed_VelocityField(veloc_Profile, param);

t_range = [t_0, t_1];
if strcmp(veloc_Profile, 'taylor_green')
    % Use sympletic integrator
    [x, y, ts_list] = generate_particle_trajectories_TaylorGreen(...
        t_range, Nts, sub_timestep, ...
        npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max,...
        veloc_fldStruct, kappa_fldStruct);
else
    disp('Numerical Scheme: Forward Euler');
    [x, y, ts_list] = generate_particle_trajectories(t_range, Nts, ...
        sub_timestep, ...
        npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max,...
        veloc_fldStruct, kappa_fldStruct);
end

trajdata_List(k).x = x;
trajdata_List(k).y = y;
trajdata_List(k).ts_list = ts_list;
trajdata_List(k).kappa_fldStruct = kappa_fldStruct;
trajdata_List(k).veloc_fldStruct = veloc_fldStruct;
trajdata_List(k).param = param;
trajdata_List(k).param.SamplingInterval_vis = SamplingInterval_vis;
trajdata_List(k).sub_timestep = sub_timestep;

end

%% Output
filepath = ['/data/tying/NumSDE/', veloc_Profile, '/PART_TRAJ/'];
[status,msg,msgID] = mkdir(filepath);
for k = 1:length(trajdata_List)
    param = trajdata_List(k).param;

    Pe = param.Pe;
    b_angle = param.b_angle;
    b_mag = param.b_mag;
    Npart = param.Npart;
    TempIntv = param.TempIntv;
    TempIntv_vis = param.TempIntv_vis;
    SamplingInterval_vis = param.SamplingInterval_vis;
    nReal = param.nReal;
    
    if (strcmp(veloc_Profile, 'taylor_green_noisy')) || (strcmp(veloc_Profile, 'childress_soward'))
        traj_filename = ['traj_', veloc_Profile, '_', kappa_Profile, ...
            '_Pe', num2str(Pe), '_Npart', num2str(Npart), ...
            '_TempIntv',num2str(TempIntv_vis), '_h', num2str(SamplingInterval_vis), ...
            '_nReal',num2str(nReal),'.mat'];
        
    elseif (strcmp(veloc_Profile, 'tg_w_mean'))
        traj_filename = ['traj_', veloc_Profile,  ...
            '_bmag', num2str(b_mag), '_bangle', num2str(b_angle), ...
            '_', kappa_Profile, '_Pe', num2str(Pe), ...
            '_Npart', num2str(Npart), '_TempIntv', num2str(TempIntv_vis), '_h', ...
            num2str(SamplingInterval_vis), '_nReal',num2str(nReal), '.mat'];
    elseif (strcmp(veloc_Profile, 'cosine'))
        traj_filename = ['traj_', veloc_Profile, ...
            '_', kappa_Profile, '_Pe', num2str(Pe), ...
            '_Npart', num2str(Npart), '_TempIntv',num2str(TempIntv_vis), '_h', ...
            num2str(SamplingInterval_vis), '_nReal',num2str(nReal), '.mat'];    
    elseif (contains(veloc_Profile, 'linear_jet'))
        traj_filename = ['traj_', veloc_Profile, ...
            '_', kappa_Profile, '_Pe', num2str(Pe), ...
            '_Npart', num2str(Npart), '_TempIntv',num2str(TempIntv_vis), '_h', ...
            num2str(SamplingInterval_vis), '_nReal',num2str(nReal), '.mat'];   
    end

    traj_fullpath = [filepath, traj_filename]
    
    % Unfolding variables
    x = trajdata_List(k).x;
    y = trajdata_List(k).y;
    ts_list = trajdata_List(k).ts_list;
    kappa_fldStruct = trajdata_List(k).kappa_fldStruct;
    veloc_fldStruct = trajdata_List(k).veloc_fldStruct;
    param = trajdata_List(k).param;
    sub_timestep = trajdata_List(k).sub_timestep;
    
    %disp('Not saving the trajectories!')
    disp('Saving the trajectories!')
    save(traj_fullpath, 'x', 'y', 'ts_list', 'kappa_fldStruct', 'veloc_fldStruct', 'param', 'sub_timestep')
end



%% Prelim Check
for  k = 1:4:length(trajdata_List)
    figure(k);
    subplot(2,2,1)
    plot(trajdata_List(k).x,trajdata_List(k).y)

    subplot(2,2,2)
    scatter(trajdata_List(k).x(end,:),trajdata_List(k).y(end,:))

    subplot(2,2,3)
    scatter(trajdata_List(k).x(1,:),trajdata_List(k).y(1,:))

    subplot(2,2,4)
    h1 = histogram(trajdata_List(k).x(end,:));
    h1.Normalization = 'probability';
    hold on
    h2 = histogram(trajdata_List(k).y(end,:));
    h2.Normalization = 'probability';
    
    
    pause
end

% Quick check for convergence
SIv_rate_vector = 1:1:4096*2;
KXX_Freq = SIv_rate_vector*0;

for i = 1:length(SIv_rate_vector)
SIv_rate = SIv_rate_vector(i);
xds = trajdata_List(k).x(1:SIv_rate:end,:);
yds = trajdata_List(k).y(1:SIv_rate:end,:);

SIv = trajdata_List(k).ts_list(1+SIv_rate) - trajdata_List(k).ts_list(1);
diffxds = diff(xds, 1);
diffyds = diff(yds, 1);
KXX_Freq(i) = var(diffxds(:))/(2*SIv);
KYY_Freq(i) = var(diffyds(:))/(2*SIv);
end

plot(SIv_rate_vector, KXX_Freq, SIv_rate_vector, KYY_Freq)

% 
% % Test Case: Shear Flow
% % Without considering the mean flow
% % i.e. Absolute diffusivity
% % Compute variance at t
% % XX, YY, XY
% cov_xy = zeros(length(ts_list), 3);
% 
% for i_t = 1:length(ts_list)
%     x_i_t = x(i_t, :);
%     y_i_t = y(i_t, :);
% 
%     cov_temp = cov(x_i_t, y_i_t);
%     
%     cov_xy(i_t, 1) = cov_temp(1,1);
%     cov_xy(i_t, 2) = cov_temp(1,2);
%     cov_xy(i_t, 3) = cov_temp(2,2);
% end
% plot(ts_list, cov_xy);
% title('Absoluate Variance of All Particles')
% 
% 
% if strcmp(veloc_Profile, 'shear')
%     A_theo = [0, GradU_scale(1); 0, 0];
%     kappa_theo = kappa_fldStruct.kappa(0, 0);
% 
%     
%     % Data: Remove trend
%     % Varying Sampling Interval to observe the transition density
%     % Compute
%     sInt_step_max = 512*2;
%     
%     sInt_list = zeros(sInt_step_max, 1);
%     Sigma_t_ts = zeros(sInt_step_max, 3);
%     
%     for sInt_step = 1:sInt_step_max
%         filter_list = 1:sInt_step:length(ts_list);
%         
%         ts_list_filtered = ts_list(filter_list);
%         x_fi = x(filter_list, :);
%         y_fi = y(filter_list, :);
%         
%         sInt = ts_list_filtered(2) - ts_list_filtered(1);
%         
%         diffx = diff(x_fi, 1)-y_fi(1:end-1,:)*GradU_scale(1)*sInt;
%         diffy = diff(y_fi, 1);
%         
%         Sigma_t = cov(diffx(:), diffy(:)); % if /(2*sInt) furthur => EA diffusivity;
%         
%         sInt_list(sInt_step) = sInt;
%         Sigma_t_ts(sInt_step, :) = [Sigma_t(1,1), Sigma_t(2,2), Sigma_t(1,2)];
%     end
%     
%     plot(sInt_list, Sigma_t_ts)
%     legend('K11', 'K22', 'K12')
%     
%     hold on
% %    % Exact kappa: divide Sigma_t_ts by (2*sIntat). Send sampling interval -> 0    
% %    plot(sInt_list, [kappa_theo(1,1), kappa_theo(2,2), kappa_theo(1,2)].*ones(length(sInt_list),1), ':' )
%   
%     
%     A = A_theo;
%     kappa = kappa_theo;
%     
%     
%     q = sqrt(-det(A));
%     Sigma_t_ts_theo = zeros(length(sInt_list), 3);
%     for i_t = 1:length(sInt_list)
%         t = sInt_list(i_t);
%         
%         if det(A) ~= 0
%             Sigma_t_theo = (t+sinh(2*q*t)/(2*q))*kappa + sinh(q*t)*sinh(q*t)/(q*q)*(A*kappa+kappa*A') + (1/q^2)*(sinh(2*q*t)/(2*q)-t)*A*kappa*A';
%         else
%             %Lik_Covar_theo = 2*kappa*[t+1/3*alpha*alpha*t^3, 0.5*alpha*t^2; 0.5*alpha*t^2, t];
%             Sigma_t_theo = 2*t*kappa + t^2*(A*kappa+kappa*A') + ((2/3)*t^3)* A*kappa*A';
%         end
%         
%         Sigma_t_ts_theo(i_t, 1) = Sigma_t_theo(1,1);
%         Sigma_t_ts_theo(i_t, 2) = Sigma_t_theo(2,2);
%         Sigma_t_ts_theo(i_t, 3) = Sigma_t_theo(1,2);
%     end
%     
%     plot(sInt_list, Sigma_t_ts_theo, ':');
%     legend('K11_T', 'K22_T', 'K12_T')
%     title('If solid line coincides with the dotted line, it means the SDE solver is working properly.')
% end