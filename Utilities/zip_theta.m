function theta_zipped  = zip_theta(theta_store)
    [Ncells, NVars, MCMC_NSample] = size(theta_store);
    
    Nsteps_pc = MCMC_NSample-1;
    
    assert(mod(Nsteps_pc, NVars) == 0);
    NSamples_zipped = Nsteps_pc/NVars;
    
    theta_zipped = zeros(Ncells, NVars, NSamples_zipped+1);
    
    % Store initial conditions
    theta_zipped(:,:, 1) = theta_store(:,:, 1);
    
    for cmp = 1:NVars
        theta_zipped(:, cmp, 2:end) = theta_store(:, cmp, (1+cmp):NVars:end);
    end
end

% % Test function
% NVars = 4;
% Ncells = 1;
% MCMC_NSample = NVars*3;
% theta_store = zeros(Ncells, NVars, MCMC_NSample + 1);
% 
% vec = theta_store(:,:, 1);
% for samp_ind = 1:(size(theta_store, 3)-1)
%     cmp = mod(samp_ind-1, NVars) + 1;
%     vec(:,cmp) = vec(:,cmp) + rand(1);
%     
%     theta_store(:,:,samp_ind+1) = vec;
% end
% 
% theta_store_vis = squeeze(theta_store)
% 
% theta_zipped  = zip_theta(theta_store);
% 
% 
% theta_zipped_vis = squeeze(theta_zipped)
