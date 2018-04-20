%% Generation of particle trajectories
%% Specification 
% Input required:

% t_0, t_1: Initial and final time
% Nts: Number of time steps between t_0 and t_1 (should be very high)
% npart_pc_DIR: refer to qg_particles_parameters.F90
% ncell_DIR:
% PART_ranges_min:
% PART_ranges_max: 

% Output:
% Python-friendly: for 2D-arary x[t, part] and y[t, part]
% x(t, part): 2D-arary; x-coordinate of part at time t
% y(t, part): 2D-arary; y-coordinate of part at time t
% ts_list: 1D-arrary (column vector); time t

%% Computation
function [x, y, ts_list] = generate_particle_trajectories_TaylorGreen(t_range, Nts, sub_timestep, npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max, veloc_fldStruct, kappa_fldStruct)
t_0 = t_range(1);
t_1 = t_range(2);

assert(t_0 < t_1);

ts_list = linspace(t_0, t_1, Nts+1)';


% Specify the SDE used
Ito = 0;
AntiIto = 1;

u_scale = veloc_fldStruct.u_scale;
kappa = kappa_fldStruct.kappa;

% Evenly distribute particles in a sub-domain
% Initialise the particle positions into two arrays x0, y0
[x_old, y_old, nparticles] = initialise_particle_positions(npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max);

x = zeros(Nts+1, nparticles);
y = zeros(Nts+1, nparticles);

x(1,:) = x_old; 
y(1,:) = y_old;

% Ref: Paliotis et al. 2009: Calculating Effective Diffusivity
e_1 = [1; 1]; e_2 = [1; -1];
d_1 = -0.5*e_2; d_2 = -0.5*e_1;

% Determine the two split velocity field
v1 = @(x, y) u_scale*sin(x);
v2 = @(x, y) u_scale*sin(x);

x_vec = [x_old; y_old];

for ts_ind = 1:Nts
    h = ts_list(ts_ind+1)-ts_list(ts_ind);
    
    % 16 sub-time step for diffusion
    %sub_timestep = 1;
    h_eff = h/sub_timestep;

    % Integrator: effectively (probably not a very good idea): 
    % x_new = S(A2(A1(x_old))))
    for h_i = 1:sub_timestep
        for part = 1:nparticles
            % A1
            x_vec(:,part) = x_vec(:,part) + h_eff*d_1*v1(dot(e_1, x_vec(:,part)));
            % A2
            x_vec(:,part) = x_vec(:,part) + h_eff*d_2*v2(dot(e_2, x_vec(:,part)));
            
%             % S
%             kappa_const = kappa(x_vec(1,part), x_vec(2,part));
%             mu = [0, 0];
%             K_sigma = 2*h_eff*kappa_const*[1, 1.5; 1.5, 3];   % Ad-hoc
%             
%             x_vec(:,part) = x_vec(:,part) + mvnrnd(mu,K_sigma,1)';            
        end
        
        % Add the noise after the deterministic advection
        % S        
        stochastic_drift = sqrt(h_eff)*sqrt(2*kappa(x_vec(1,:), x_vec(2,:))).*randn([2, nparticles]);
        
        x_vec = x_vec + stochastic_drift;
    end
        
    x(ts_ind+1,:) = x_vec(1,:);
    y(ts_ind+1,:) = x_vec(2,:);
end

end