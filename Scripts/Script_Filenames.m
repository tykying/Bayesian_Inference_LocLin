DataFolder_BI = '/data/tying/BayesianInference/';
DataFolder_Video = '/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Video_Generated/';
DataFolder_Figure = '/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Figure_Generated/';
DataFolder_Publication = '/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Figure_Publications/';


if contains(veloc_Profile, 'QGM2')
    veloc_Profile_Folder = ['/Layer_', num2str(Layer), '/' veloc_Profile, '/'];
else
    veloc_Profile_Folder = ['/' veloc_Profile, '/'];
end

if strcmp(DataTypeString, 'BinnedTrajData')
    % Dependent on: 
    % 1) veloc_Profile 2)DS_rate 3) SamplingInterval_vis 4) Nx_cell_ARG
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_BTD = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'BinnedTrajTransMat')
    % Dependent on: (SAME AS 'BinnedTrajData')
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_BTTM = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'BayesSampleStat')
    % Dependent on:
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    % 5) InfScheme  6) Nsteps_pcv
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    DataParam = [DataParam, '_', InfScheme, '_Nsteps', num2str(Nsteps_pcv)];
    
    % Ad-hoc for tg+mean drift
    if strcmp(veloc_Profile, 'tg_w_mean')
        DataParam = [DataParam, '_uscale', num2str(u_scale), '_uangle', num2str(u_angle)];
    end
    
    DataParam = [DataParam, '_nTrial', num2str(nTrial)];

    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_BSS = [filepath, DataTypeString, '_', DataParam, '.mat'];
    
    filename_BSSPostSpin = [filepath, DataTypeString, 'PostSpin_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'BayesMHRW')
    % Dependent on:
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    % 5) InfScheme  6) Nsteps_pcv
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    DataParam = [DataParam, '_', InfScheme];
    
    DataParam = [DataParam, '_nTrial', num2str(nTrial)];

    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_BMHRW = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'BayesMHConfig')
    % Dependent on:
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    % 5) InfScheme  6) Nsteps_pcv
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    DataParam = [DataParam, '_', InfScheme];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_BMConfig = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'BayesSampleMeanFine')
    % Dependent on: (SAME AS 'BayesSampleStat')
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    % 5) InfScheme  6) Nsteps_pcv
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    DataParam = [DataParam, '_', InfScheme, '_Nsteps', num2str(Nsteps_pcv)];
        
    % Ad-hoc for tg+mean drift
    if strcmp(veloc_Profile, 'tg_w_mean')
        DataParam = [DataParam, '_uscale', num2str(u_scale), '_uangle', num2str(u_angle)];
    end
    
    DataParam = [DataParam, '_nTrial', num2str(nTrial)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_BSMF = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'ParticleBinnedEvol')
    % Dependent on: (SAME AS 'BayesSampleStat')
    % 1) Layer  2) SamplingInterval_vis % 3) Nx_cell_ARG 
    DataParam = [veloc_Profile, '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];    
    
    filepath = [DataFolder_Video, DataTypeString, veloc_Profile_Folder];
    filename_PBE = [filepath, DataTypeString, '_', DataParam, '.avi'];
end


if strcmp(DataTypeString, 'MeanSampIntvConv')
    % Dependent on: (SAME AS 'BayesSampleStat')
    % 1) Layer  2) SamplingInterval_vis % 3) Nx_cell_ARG 
    DataParam = [FieldName, '_', veloc_Profile, '_Idx', num2str(Nx_cell_ARG)];    
    
    filepath = [DataFolder_Figure, DataTypeString, veloc_Profile_Folder];
    filename_MSIC = [filepath, DataTypeString, '_', DataParam, '.', ImageFormat(1:3)];
end


if strcmp(DataTypeString, 'FlowDecomp')
    % Dependent on: (SAME AS 'BinnedTrajData')
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_FD = [filepath, DataTypeString, '_', DataParam, '.mat'];
end


if strcmp(DataTypeString, 'ACTransStat')
    % Dependent on: (SAME AS 'BinnedTrajData')
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_ACTS = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'ACCorreStat')
    % Dependent on: (SAME AS 'BinnedTrajData')
    % 1) veloc_Profile  2)DS_rate  3) SamplingInterval_vis % 4) Nx_cell_ARG 
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_ACCS = [filepath, DataTypeString, '_', DataParam, '.mat'];
end

if strcmp(DataTypeString, 'RelDiffSampIntConv')
    DataParam = [FieldName, '_Idx', num2str(Nx_cell_ARG)];    
    
    filepath = [DataFolder_Figure, DataTypeString, veloc_Profile_Folder];
    filename_RDSIC = [filepath, DataTypeString, '_', DataParam, '.', ImageFormat(1:3)];
end


if strcmp(DataTypeString, 'KvsUData')
    DataParam = [veloc_Profile, '_DSr', num2str(DS_rate), '_h', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell_ARG)];
    
    filepath = [DataFolder_BI, DataTypeString, veloc_Profile_Folder];
    filename_ACCS = [filepath, DataTypeString, '_', DataParam, '.mat'];
end