%% Load Trajectory Data Script
if (~contains(veloc_Profile, 'QGM2_DS')==0)  % =0: Matched
    % Old Method:
    % Large amount of data setting (move to _HD file)
    Npart = 2^20;   % Hardcoded
    TempIntv_vis = 3653;  % Hardcoded
    
    % Ad-hoc: Find out the right file name;
    % 8: starting from h=8, all particle trajectories are stored
    SamplingInterval_vis_file = min(SamplingInterval_vis, 8);
    Npart_file = min([Npart/(8/SamplingInterval_vis_file), Npart]);   % Hardcoded
    TempIntv_vis_file = TempIntv_vis - mod(TempIntv_vis, SamplingInterval_vis_file);  % Hardcoded
    
    % New Method:
    % 30 by 30 realistic setting
    Npart_file = 30*30;   % Hardcoded
    TempIntv_vis_file = 3653;  % Hardcoded
    SamplingInterval_vis_file = 1;  % Hardcoded
        
    filepath = ['/data/tying/QGM2/PART_Layer_', num2str(Layer), '/PART_TRAJ/'];
    traj_filename = ['traj_Npart', num2str(Npart_file), ...
                    '_TempIntv', num2str(TempIntv_vis_file), ...
                    '_h',num2str(SamplingInterval_vis_file), '_nReal1.mat'];

    traj_fullpath = [filepath, traj_filename]
    
    if strcmp(traj_fullpath_old, traj_fullpath)
        disp('No need to reload the Trajectories.')
    else
        load(traj_fullpath);
        
        % Read x,y,ts_list,L,ngrid
        %%% Convert x, y from centimeters to meters
        x = x/(100); y = y/(100);  L = L/(100);   % in meters
        h = (ts_list(2)-ts_list(1));  % in seconds
        
%         time_scale = 3600*24;
%         space_scale = 1000;
        time_scale = 1;
        space_scale = 1;
        
        traj_fullpath_old = traj_fullpath;
    end
end

if (~contains(veloc_Profile, 'QGM2_nReal_DS')==0)
    Npart = 1048576;   % = 2^20; Hardcoded
    nReal = 60;   % Hardcoded
    TempIntv_vis = 4*SamplingInterval_vis;  % Hardcoded

    filepath = ['/data/tying/QGM2/PART_Layer_', num2str(Layer), '/PART_TRAJ/'];
    traj_filename = ['traj_Npart', num2str(Npart), '_TempIntv',num2str(TempIntv_vis), '_h',num2str(SamplingInterval_vis), '_nReal',num2str(nReal),'.mat'];
    traj_fullpath = [filepath, traj_filename]
    
    if strcmp(traj_fullpath_old, traj_fullpath)
        disp('No need to reload the Trajectories.')
    else
        load(traj_fullpath);
        
        % Read x,y,ts_list,L,ngrid
        %%% Convert x, y from centimeters to meters
        x = x/(100); y = y/(100);  L = L/(100);   % in meters
        h = (ts_list(2)-ts_list(1));  % in seconds
        
        time_scale = 1;
        space_scale = 1;
        
        traj_fullpath_old = traj_fullpath;
    end
end


if (strcmp(veloc_Profile, 'shear') || strcmp(veloc_Profile, 'linear_shear'))
    kappa_Profile = 'const_skew';
    nReal = 1;
    SamplingInterval_vis_file = 1;  % Hardcoded
    Npart = 2^6;
    TempIntv_vis = 1;

    
    filepath = ['/data/tying/NumSDE/', veloc_Profile, '/PART_TRAJ/'];
    traj_filename = ['traj_', veloc_Profile, '_', kappa_Profile, ...
                    '_Pe', num2str(Pe), '_Npart', num2str(Npart), ...
                    '_TempIntv',num2str(TempIntv_vis), '_h', num2str(SamplingInterval_vis_file), ...
                    '_nReal',num2str(nReal),'.mat'];
    traj_fullpath = [filepath, traj_filename]

    load(traj_fullpath);
    
    % Read x,y,ts_list,L,ngrid
    %%% Convert x, y from centimeters to meters
    h = (ts_list(2)-ts_list(1));  % in seconds
    
    time_scale = 1;
    space_scale = 1;
end

if contains(veloc_Profile, 'linear_jet')
    kappa_Profile = 'const_skew';
    nReal = 1;
    SamplingInterval_vis_file = 0;  % Hardcoded
    Npart = 2^8;
    %Npart = 2^10;
    TempIntv_vis = 32;

    
    filepath = ['/data/tying/NumSDE/', veloc_Profile, '/PART_TRAJ/'];
    traj_filename = ['traj_', veloc_Profile, '_', kappa_Profile, ...
                    '_Pe', num2str(Pe), '_Npart', num2str(Npart), ...
                    '_TempIntv',num2str(TempIntv_vis), '_h', num2str(SamplingInterval_vis_file), ...
                    '_nReal',num2str(nReal),'.mat'];
    traj_fullpath = [filepath, traj_filename]

    load(traj_fullpath);
    
    % Read x,y,ts_list,L,ngrid
%     %%% Convert x, y from centimeters to meters
%     x = x*1E-5; y = y*1E-5;
%     ts_list = ts_list/(24*3600);
%     
%     h = (ts_list(2)-ts_list(1));  % in days
%     
%     time_scale = 3600*24;
%     space_scale = 1000;    
    h = (ts_list(2)-ts_list(1));  % in seconds
    
    time_scale = 1;
    space_scale = 1;
end



if (strcmp(veloc_Profile, 'taylor_green_noisy')) || (strcmp(veloc_Profile, 'childress_soward'))
    %kappa_Profile = 'const';   % Hardcoded
    nReal = 1;   % Hardcoded
    SamplingInterval_vis_file = 1;  % Hardcoded
    %Pe = 100;  % Hardcoded
    %disp('Hardcode Pe=100 in LoadTrajData');
    if (strcmp(veloc_Profile, 'taylor_green_noisy'))
        Npart = 2^16; % Hardcoded
    elseif (strcmp(veloc_Profile, 'childress_soward'))
        Npart = 2^16;
    end
    
    TempIntv_vis = 16; % Hardcoded

    
    filepath = ['/data/tying/NumSDE/', veloc_Profile, '/PART_TRAJ/'];
    traj_filename = ['traj_', veloc_Profile, '_', kappa_Profile, ...
                    '_Pe', num2str(Pe), '_Npart', num2str(Npart), ...
                    '_TempIntv',num2str(TempIntv_vis), '_h', num2str(SamplingInterval_vis_file), ...
                    '_nReal',num2str(nReal),'.mat'];
    traj_fullpath = [filepath, traj_filename]

    load(traj_fullpath);
    
    % Read x,y,ts_list,L,ngrid
    %%% Convert x, y from centimeters to meters
    h = (ts_list(2)-ts_list(1));  % in seconds
    
    time_scale = 1;
    space_scale = 1;
end


if (strcmp(veloc_Profile, 'tg_w_mean')) 
    %kappa_Profile = 'const';   % Hardcoded
    nReal = 1;   % Hardcoded
    SamplingInterval_vis_file = 0.03125;  % Hardcoded   
    
    if contains(RunProfile, '_TEST') == 1
        Npart = 2^8; % Hardcoded
    else
        Npart = 2^10; % Hardcoded
    end
    
    TempIntv_vis = 256; % Hardcoded

    
    filepath = ['/data/tying/NumSDE/', veloc_Profile, '/PART_TRAJ/'];
        traj_filename = ['traj_', veloc_Profile, ...
            '_uscale', num2str(u_scale), '_uangle', num2str(u_angle), ...
            '_', kappa_Profile, '_Pe', num2str(Pe), ...
            '_Npart', num2str(Npart), '_TempIntv',num2str(TempIntv_vis), '_h', ...
            num2str(SamplingInterval_vis_file), '_nReal',num2str(nReal) '.mat'];
        
    traj_fullpath = [filepath, traj_filename]

    load(traj_fullpath);
    
    % Read x,y,ts_list,L,ngrid
    %%% Convert x, y from centimeters to meters
    h = (ts_list(2)-ts_list(1));  % in seconds
    
    time_scale = 1;
    space_scale = 1;
end


if (strcmp(veloc_Profile, 'cosine')) 
    %kappa_Profile = 'const';   % Hardcoded
    nReal = 1;   % Hardcoded
    
    TempIntv_vis = 64; % Hardcoded

    SamplingInterval_vis_file = TempIntv_vis/2^13;  % Hardcoded   
    
    Npart = 2^12; % Hardcoded
    
    
    filepath = ['/data/tying/NumSDE/', veloc_Profile, '/PART_TRAJ/'];
        traj_filename = ['traj_', veloc_Profile, '_', kappa_Profile, '_Pe', num2str(Pe), ...
            '_Npart', num2str(Npart), '_TempIntv',num2str(TempIntv_vis), '_h', ...
            num2str(SamplingInterval_vis_file), '_nReal',num2str(nReal), '.mat'];
        
    traj_fullpath = [filepath, traj_filename]

    load(traj_fullpath);
    
    % Read x,y,ts_list,L,ngrid
    %%% Convert x, y from centimeters to meters
    h = (ts_list(2)-ts_list(1));  % in seconds
    
    time_scale = 1;
    space_scale = 1;
end


if strcmp(veloc_Profile, 'Mixed')
    % Mix two types of trajectories together
    veloc_Profile = 'childress_soward';
    kappa_Profile = 'sinusoidal';
    run('Script_LoadTrajData');
    x_cs = x; y_cs = y; ts_list_cs = ts_list;
    
    veloc_Profile = 'taylor_green_noisy';
    kappa_Profile = 'const';
    run('Script_LoadTrajData');
    x_tg = x; y_tg = y; ts_list_tg = ts_list;
    
    x = [x_cs, x_tg];
    y = [y_cs, y_tg];
    assert(all(size(ts_list_cs) == size(ts_list_tg)));
    
    veloc_Profile = 'Mixed';
    disp('Finished Loading Mixed Trajectory Data.')
end

