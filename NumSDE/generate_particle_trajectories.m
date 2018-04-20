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
function [x, y, ts_list] = generate_particle_trajectories(t_range, Nts, sub_timestep, npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max, veloc_fldStruct, kappa_fldStruct)
t_0 = t_range(1);
t_1 = t_range(2);

assert(t_0 < t_1);

ts_list = linspace(t_0, t_1, Nts+1)';


% Specify the SDE used
Ito = 0;
AntiIto = 1;

% Evenly distribute particles in a sub-domain
% Initialise the particle positions into two arrays x0, y0
[x_old, y_old, nparticles] = initialise_particle_positions(npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max);

x = zeros(Nts+1, nparticles);
y = zeros(Nts+1, nparticles);

x(1,:) = x_old; y(1,:) = y_old;

for ts_ind = 1:Nts
    h = ts_list(ts_ind+1)-ts_list(ts_ind);
    
    [x_new, y_new] = timestepping_FEuler(x_old,y_old, h, sub_timestep, veloc_fldStruct, kappa_fldStruct, AntiIto);

    %[x_new, y_new] = impose_reflectiveBC(x_new, y_new, PART_ranges_min, PART_ranges_max);
    
    x(ts_ind+1, :) = x_new;
    y(ts_ind+1, :) = y_new;
    
    x_old = x_new; y_old = y_new;
    
    if mod(log2(ts_ind), 2)==1
        disp(['SDE time step = ', num2str(ts_ind)])
    end

end

end