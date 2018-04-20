%% Initialise the struct
% Specifications:

% Input:
% x, y: trajectories of all particles
% ts_list: time step

% Output:
% Jumps(:): struct array that contains all individual jumps
%% Computation

function [Jumps, TrajData_SubSampled] = convert_traj_to_JumpsStruct(x, y, ts_list, SamplingInterval, SamplingParticle)

TrajData_SubSampled = SubSample_Traj(x, y, ts_list, SamplingInterval, SamplingParticle);

Jumps = struct('x_t0',{}, 'x_t1', {}, 'h', {}, 'dx', {});

x = TrajData_SubSampled.x_SubSampled;
y = TrajData_SubSampled.y_SubSampled;
ts_list = TrajData_SubSampled.ts_list_SubSampled;

diffx = diff(x,1);
diffy = diff(y,1);

h = diff(TrajData_SubSampled.ts_list_SubSampled);

% Make sure there are at least 2 samples
assert(length(ts_list) >= 2);
assert(length(h) == size(diffx, 1));

% Assign the struct
jump = 0;
for t_ind = 1:length(h)
    for part = 1:size(x, 2)
        jump = jump + 1;
        
        Jumps(jump).x_t0 = [x(t_ind, part); y(t_ind, part)];
        Jumps(jump).x_t1 = [x(t_ind+1, part); y(t_ind+1, part)];
        Jumps(jump).dx   = [diffx(t_ind, part); diffy(t_ind, part)];

        Jumps(jump).h    = h(t_ind);
    end
end

end