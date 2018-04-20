%% Setting-up Parameters for Down-Sampling
DS_frac = 1/DS_rate;


%% For Linear SDE
if ( strcmp(veloc_Profile, 'shear') ...
        || strcmp(veloc_Profile, 'linear_shear') ...
        || strcmp(veloc_Profile, 'taylor_green_noisy') ...
        || strcmp(veloc_Profile, 'childress_soward')  ...
        || strcmp(veloc_Profile, 'Mixed') )

    % IDEA: 
    % For NJumps_target=2^(2n), varies number of particles and sampling interval
    % such that NJumps = 2^n, n = max_level
    
    % Ensure NJumps_RAW = 2^(2n)
    NJumps_target = NJumps_RAW*DS_frac;
    %assert(mod(log2(NJumps_target), 2) == 0)
    
    max_level = log2(NJumps_target)/2;  % n = max_level in the IDEA
    level = log2(SamplingInterval_vis);
    
    SamplingInterval = 2^(level)*h;
    SamplingParticle = 2^(max_level-level);
    
    % Ad-hoc: fix to use all particles
    SamplingInterval = 2^(level)*h;
    SamplingParticle = DS_rate;
    
    
    t_offset = ts_list(1);  % Start at the beginning
    
    % Same as the Naive QGM2
    % Downsample ts_list
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer
    
    % Downsample particles
    SampPart_DInd = SamplingParticle;
    
    % Output Indices for Down-Sampling
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    DS_t_Ind = t_Ind_offset:SampIntv_DInd:Nts_RAW;
    DS_Part_Ind = 1:SampPart_DInd:nparticles_RAW;
    
    TrajData = struct('x', x(DS_t_Ind, DS_Part_Ind), ...
                      'y', y(DS_t_Ind, DS_Part_Ind), ...
                      'ts_list', ts_list(DS_t_Ind));
    
    % N.B. Also can use `subsample_Trajectories', which does not assume uniform time step
    % Need to filter out x, y, ts_list for offset though
    % TrajData = subsample_Trajectories(x, y, ts_list, SamplingInterval, SamplingParticle);
end

if ( strcmp(veloc_Profile, 'tg_w_mean'))
    % IDEA: 
    % DS_rate: Use all particles; but vary the sampling interval
    
    SamplingInterval = SamplingInterval_vis*h;
    SamplingParticle = DS_rate;    
    
    t_offset = ts_list(1);  % Start at the beginning
    
    % Same as the Naive QGM2
    % Downsample ts_list
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer
    
    % Downsample particles
    SampPart_DInd = SamplingParticle;
    
    % Output Indices for Down-Sampling
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    DS_t_Ind = t_Ind_offset:SampIntv_DInd:Nts_RAW;
    DS_Part_Ind = 1:SampPart_DInd:nparticles_RAW;
    
    TrajData = struct('x', x(DS_t_Ind, DS_Part_Ind), ...
                      'y', y(DS_t_Ind, DS_Part_Ind), ...
                      'ts_list', ts_list(DS_t_Ind));
    
    % N.B. Also can use `subsample_Trajectories', which does not assume uniform time step
    % Need to filter out x, y, ts_list for offset though
    % TrajData = subsample_Trajectories(x, y, ts_list, SamplingInterval, SamplingParticle);
end


if ( strcmp(veloc_Profile, 'cosine'))
    % IDEA: 
    % DS_rate: Use all particles; but vary the sampling interval
    
    SamplingInterval = SamplingInterval_vis*h;
    SamplingParticle = DS_rate;    
    
    t_offset = ts_list(1);  % Start at the beginning
    
    % Same as the Naive QGM2
    % Downsample ts_list
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer
    
    % Downsample particles
    SampPart_DInd = SamplingParticle;
    
    % Output Indices for Down-Sampling
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    DS_t_Ind = t_Ind_offset:SampIntv_DInd:Nts_RAW;
    DS_Part_Ind = 1:SampPart_DInd:nparticles_RAW;
    
    TrajData = struct('x', x(DS_t_Ind, DS_Part_Ind), ...
                      'y', y(DS_t_Ind, DS_Part_Ind), ...
                      'ts_list', ts_list(DS_t_Ind));
    
    % N.B. Also can use `subsample_Trajectories', which does not assume uniform time step
    % Need to filter out x, y, ts_list for offset though
    % TrajData = subsample_Trajectories(x, y, ts_list, SamplingInterval, SamplingParticle);
end

if ( contains(veloc_Profile, 'linear_jet') )
    % IDEA: 
    % DS_rate: Use all particles; but vary the sampling interval
    
    SamplingInterval = SamplingInterval_vis*24*3600/time_scale;
    SamplingParticle = DS_rate;    
    
    t_offset = ts_list(1);  % Start at the beginning
    
    % Same as the Naive QGM2
    % Downsample ts_list
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer
    
    % Downsample particles
    SampPart_DInd = SamplingParticle;
    
    % Output Indices for Down-Sampling
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    DS_t_Ind = t_Ind_offset:SampIntv_DInd:Nts_RAW;
    DS_Part_Ind = 1:SampPart_DInd:nparticles_RAW;
    
    TrajData = struct('x', x(DS_t_Ind, DS_Part_Ind), ...
                      'y', y(DS_t_Ind, DS_Part_Ind), ...
                      'ts_list', ts_list(DS_t_Ind));
    
    % N.B. Also can use `subsample_Trajectories', which does not assume uniform time step
    % Need to filter out x, y, ts_list for offset though
    % TrajData = subsample_Trajectories(x, y, ts_list, SamplingInterval, SamplingParticle);
end


%% For QGM2 (Long Trajectories)
if strcmp(veloc_Profile, 'QGM2_DSpart') || strcmp(veloc_Profile, 'QGM2_DStemp')
    % IDEA: DS_frac controls number of particles
    % Naive DownSampling: 
    % Given Sampling Interval and DS_frac,
    % Down-sample ONLY the number of particles using DS_frac
    % Problem: NJumps decreases as Sampling Interval increases
    TempIntv = (ts_list(end) - ts_list(1));
    SamplingInterval = 24*3600*SamplingInterval_vis;
    
    % t_offset = ts_list(1) + 6 * (365*24*3600);  % Start at the 6th year
    t_offset = ts_list(1);  % Start at the beginning
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    
    % Downsample ts_list; DInd = Delta Index
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer
    
    % Downsample particles (Naive; without preserving total NJumps)
    SampPart_DInd = round(1/DS_frac);  % Should be an integer

    % If DS_rate > 1: 
    % preserve the total number of particles at short sampling interval h
    if (DS_rate > 1) && (SamplingInterval_vis <= 4)
      SampPart_DInd_Multipler = 2^(3 - log2(SamplingInterval_vis));
      SampPart_DInd = SampPart_DInd / SampPart_DInd_Multipler;
    end

   % Output Indices for Down-Sampling
    DS_t_Ind = t_Ind_offset:SampIntv_DInd:Nts_RAW;
    DS_Part_Ind = 1:SampPart_DInd:nparticles_RAW;
    
    
    TrajData = struct('x', x(DS_t_Ind, DS_Part_Ind), ...
                      'y', y(DS_t_Ind, DS_Part_Ind), ...
                      'ts_list', ts_list(DS_t_Ind));
end


if strcmp(veloc_Profile, 'QGM2_DSmSep')
    % Sketch of time line
    % mSep:    <..........>       (Set in this function) 
    % SampInt: <===>              (Given)
    % t_ind:   |---|---|---|---|---|---|---|
    % t_i:     Y---N---N---Y---N---N---Y---N
    % t_f:     N---Y---N---N---Y---N---N---Y
    % Y: Record position
    
    % DS_frac: Downsampling for Particles
    
    % IDEA: mSep: minimum separation between each jumps
    % Naive DownSampling: 
    % Given Sampling Interval, DS_frac and mSep
    % Down-sample ONLY the number of particles using DS_frac
    % Problem: NJumps decreases as Sampling Interval increases
    TempIntv = (ts_list(end) - ts_list(1));
    SamplingInterval = 24*3600*SamplingInterval_vis;
    
    mSep = 128*24*3600;  % mSep: minimum temporal separation between each jump
    
    % t_offset = ts_list(1) + 6 * (365*24*3600);  % Start at the 6th year
    t_offset = ts_list(1);  % Start at the beginning
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    
    % mSep_DInd: The Delta Index corresponding to separations between jumps
    mSep_DInd = round(mSep/h);
    DS_ti_Ind = t_Ind_offset:mSep_DInd:(Nts_RAW-mSep_DInd);
    
    % DInd = Delta Index
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer>=1
    DS_tf_Ind = DS_ti_Ind+SampIntv_DInd;
    
    % Downsample particles (Naive; without preserving total NJumps)
    SampPart_DInd = round(1/DS_frac); 
    DS_Part_Ind = 1:SampPart_DInd:nparticles_RAW;             

    x_i = x(DS_ti_Ind, DS_Part_Ind);
    x_f = x(DS_tf_Ind, DS_Part_Ind);
    y_i = y(DS_ti_Ind, DS_Part_Ind);
    y_f = y(DS_tf_Ind, DS_Part_Ind);
    
    % Ensure Exact Time is used
    assert(abs(SamplingInterval-(ts_list(DS_tf_Ind(1))-ts_list(DS_ti_Ind(1)))) < 1);
    
    % Merge x_i, x_f by seeing each jumps as new particles
    TrajData = struct('x', [x_i(:), x_f(:)]', ...
                      'y', [y_i(:), y_f(:)]', ...
                      'ts_list', [0, SamplingInterval]);        
end

% TODO: Fix the total number of Jumps using DS_frac
% i.e. Varies nparticles and SamplingInterval together, given a DS_frac
% -> Solve the problem of SamplingInverval increases -> Reduction of NJumps

%% For QGM2_nREAL (Multiple Realisations of Short Trajectories)
% QGM2_nREAL: ALREADY DOWNSAMPLED IN TIME
if strcmp(veloc_Profile, 'QGM2_nReal_DSSpat')
    % Sketch of particles
    % For all time t_i
    % t_1:        p1---p2---p3---p4---p5---p6---p7---p8
    % t_2:        p1---p2---p3---p4---p5---p6---p7---p8
    % Record:     YY---NN---NN---YY---NN---NN---YY---NN
    % Y: Record position
    
    % DS_frac: Downsampling for Particles
    
    % IDEA: DS_frac controls the number of spatial seed points 
    % Given Sampling Interval and DS_frac,
    % Downsample ONLY the number of particles using DS_frac
    TempIntv = (ts_list(end) - ts_list(1));
    SamplingInterval = 24*3600*SamplingInterval_vis;
    assert(abs(4*SamplingInterval-TempIntv) < 1)
    
    % t_offset = ts_list(1) + 6 * (365*24*3600);  % Start at the 6th year
    t_offset = ts_list(1);  % Start at the beginning
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);
    
    % Downsample ts_list; DInd = Delta Index
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer
    
    % ATTENTION HERE
    % Downsample particles (Need special attention to ensure the seeds points are the same at different time)
    SampPart_DInd = round(1/DS_frac);  % Should be an integer
    assert(Npart*nReal == nparticles_RAW)
    DS_Part_Ind_step = 1:SampPart_DInd:Npart;  % Seed point at each time step
    
    % Trick to collect all the seed points at all time steps
    % e.g. A = [0:(5-1)]*4 .* ones(3, 1) + [1:3]'
    DS_Part_Ind_Array = ([0:(nReal-1)]*Npart) .* ones(length(DS_Part_Ind_step), 1) + DS_Part_Ind_step';
    
    DS_Part_Ind = DS_Part_Ind_Array(:);

%     % For debugging
%     x_temp = x(:, DS_Part_Ind);
%     y_temp = y(:, DS_Part_Ind);
%     assert( all(x_temp(1, 1:10) == x_temp(1, (length(DS_Part_Ind_step)+[1:10]))) )
%     assert( all(y_temp(1, 1:10) == y_temp(1, (length(DS_Part_Ind_step)+[1:10]))) )

    TrajData = struct('x', x(:, DS_Part_Ind), ...
                      'y', y(:, DS_Part_Ind), ...
                      'ts_list', ts_list(:));
end

% QGM2_nREAL: ALREADY DOWNSAMPLED IN TIME
if strcmp(veloc_Profile, 'QGM2_nReal_DSTemp')
    % Sketch of all scenarios
    % Delta t = already equal to SamplingInterval
    % DS_frac=1 <=> Y=(t_1, t_2, ... t_5) => Each Particle has 4 jumps
    % DS_frac=0.5 <=> Y=(t_1, t_2, t_3)   => Each Particle has 2 jumps
    % DS_frac=0.25 <=> Y=(t_1, t_2)       => Each Particle has 1 jumps
    
    % N.B. When downsampling
    % length(ts_list) == 1 + 1/min(DS_frac)
    % t_Ind_end = 1 + (length(ts_list)-1/DS_frac)
    
    % IDEA: DS_frac controls the number of consecutive jumps  
    % Given Sampling Interval and DS_frac,
    % Data Loaded has been downsampled already: 
    % ts_list has already been downsampled to be length 5 (with 4 jumps)
    
    % Assumed DS_frac can only be 1, 0.5 or 0.25
    assert((DS_frac == 1) || (DS_frac==0.5) || (DS_frac==0.25));

    TempIntv = (ts_list(end) - ts_list(1));
    SamplingInterval = 24*3600*SamplingInterval_vis;
    
    % t_offset = ts_list(1) + 6 * (365*24*3600);  % Start at the 6th year
    t_offset = ts_list(1);  % Start at the beginning
    t_Ind_offset = find(ts_list >= t_offset*(1-10^-10), 1);

    if DS_frac == 0.25
        t_Ind_end = 2;
    elseif DS_frac == 0.5
        t_Ind_end = 3;
    elseif DS_frac == 1
        t_Ind_end = 5;
    end
%     t_Ind_end = 1 + (length(ts_list)-1/DS_frac);
    
    % Downsample ts_list; DInd = Delta Index
    SampIntv_DInd = round(SamplingInterval/h);  % Should be an integer=1

    % Debugging
    assert(t_Ind_offset == 1);
    assert(SampIntv_DInd == 1);
    
    DS_t_Ind = t_Ind_offset:SampIntv_DInd:t_Ind_end;
    % In fact, DS_t_Ind = 1:1:t_Ind_end

    TrajData = struct('x', x(DS_t_Ind, :), ...
                      'y', y(DS_t_Ind, :), ...
                      'ts_list', ts_list(DS_t_Ind));
end
