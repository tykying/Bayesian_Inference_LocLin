function [TrajData_SubSampled] = subsample_Trajectories(x, y, ts_list, SamplingInterval, SamplingParticle)

assert(all(size(x) == size(y)))

[NObs, nparticles] = size(x);

assert(NObs == length(ts_list));

% Assumption on Data: 
% 1) All particles share the same ts_list 
% 2) Non-uniform timestep allowed
% 3) The Sampling Interval enables at least 10 subsamples 

% Method 1: 
% % Problem of this method: some samples with with time interval in between
% are possible, e.g. ts_list = [0, 1.9, 3.0] will all be accepted with SamplingInterval = 1.5
% % Matching the timestep with the SamplingInterval
% timestep_ceil = floor((ts_list-ts_list(0))/SamplingInterval);
% Obs_list = 1:NObs;
% 
% Obs_list_sample = [true, (diff(timestep_ceil) > 0)]

% Method 2: 
% Sequentially choose the sample: 
% Interval between each subsample has minimally SamplingInterval
tn_last = ts_list(1);
Obs_list = [1];

for t_ind = 1: NObs
    tn = ts_list(t_ind);
    
    if (tn >= (tn_last + SamplingInterval)*(1-1E-10))  % Prevent machine error
        Obs_list = [Obs_list, t_ind];
        tn_last = tn;
    end
end

% SubSample Particle Trajectories
Part_list = 1:SamplingParticle:nparticles;

x_SubSampled = x(Obs_list, Part_list);
y_SubSampled = y(Obs_list, Part_list);
ts_list_SubSampled = ts_list(Obs_list)';

TrajData_SubSampled = struct('x', x_SubSampled, 'y', y_SubSampled, 'ts_list', ts_list_SubSampled);

end