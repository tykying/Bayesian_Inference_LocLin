function [theta_MAP_Sessions, logPost_MAP_Sessions] = locate_theta_MAP_Sessions(theta_store, logPost_ts, NSessions)
MCMC_NSample = size(theta_store, 3);

if NSessions > 1
    theta_MAP_Sessions = zeros( [size(theta_store, 1), size(theta_store, 2), NSessions] );
    logPost_MAP_Sessions = zeros( [size(theta_store, 1), NSessions] );
    
    for Ses = 1:NSessions
        Ses_bgn = round(MCMC_NSample*(Ses-1)/NSessions)+1;
        Ses_end = round(MCMC_NSample*Ses/NSessions);
        
        logPost_ts_Ses = logPost_ts(:, Ses_bgn:Ses_end);
        theta_store_Ses = theta_store(:, :, Ses_bgn:Ses_end);
        
        [theta_MAP, logPost_MAP, MAP_SampleInd] = locate_theta_MAP(theta_store_Ses, logPost_ts_Ses);
        
        theta_MAP_Sessions(:, :, Ses) = theta_MAP;
        logPost_MAP_Sessions(:, Ses) = logPost_MAP;
    end
    
elseif NSessions == 1
    [theta_MAP, logPost_MAP, MAP_SampleInd] = locate_theta_MAP(theta_store, logPost_ts);

    theta_MAP_Sessions = theta_MAP;
    logPost_MAP_Sessions = logPost_MAP;
end

end