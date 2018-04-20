function Homo_time_scale_ind = estimate_homo_timescale_SpatHomo(AutoCorr_SamplingInterval)
% Estimate the time for the autocorrelation with step 1 to decay to 0.05
Dn = 1;
AutoCorr_XuXs_Dn = AutoCorr_SamplingInterval(1+Dn,1,:);
AutoCorr_XuYs_Dn = AutoCorr_SamplingInterval(1+Dn,2,:);
AutoCorr_YuXs_Dn = AutoCorr_SamplingInterval(1+Dn,3,:);
AutoCorr_YuYs_Dn = AutoCorr_SamplingInterval(1+Dn,4,:);
AutoCorr_Dn1 = [AutoCorr_XuXs_Dn(:), AutoCorr_XuYs_Dn(:), AutoCorr_YuXs_Dn(:), AutoCorr_YuYs_Dn(:)];

% 3b. Estimate the time for the autocorrelation with step 1 to decay to 0.05
% Find from the back the first false: AutoCorr_Dn1(:,1) < 0.05
threshold_ind_RuRs = zeros(4, 1);
for i = 1:4    % 4: _XuXs, _XuYs, _YuXs, _YuYs
    threshold_ind_RuRs_i = find((AutoCorr_Dn1(:,i) > 0.05), 1, 'last');
    if isempty(threshold_ind_RuRs_i)
        threshold_ind_RuRs(i) = 1;
    else
        threshold_ind_RuRs(i) = threshold_ind_RuRs_i;
    end
end
threshold_ind = max(threshold_ind_RuRs);


% Ensure that when step Dn = 2: the Autocorrelation has also decayed
Dn = 2;
AutoCorr_XuXs_Dn = AutoCorr_SamplingInterval(1+Dn,1,:);
AutoCorr_XuYs_Dn = AutoCorr_SamplingInterval(1+Dn,2,:);
AutoCorr_YuXs_Dn = AutoCorr_SamplingInterval(1+Dn,3,:);
AutoCorr_YuYs_Dn = AutoCorr_SamplingInterval(1+Dn,4,:);
AutoCorr_Dn2 = [AutoCorr_XuXs_Dn(:), AutoCorr_XuYs_Dn(:), AutoCorr_YuXs_Dn(:), AutoCorr_YuYs_Dn(:)];

assert(all(AutoCorr_Dn2(threshold_ind, :) < 0.025));


Homo_time_scale_ind = threshold_ind;
end