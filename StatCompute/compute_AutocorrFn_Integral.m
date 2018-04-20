function [AutoIntegral_lag] = compute_AutocorrFn_Integral(AutocorrFn, SamplingInterval)
[NCells, NTimeLag, NAutocorrFn] = size(AutocorrFn);

% 5 Components: Kxx, Kyy, Kxy, Kyx; Kdsds
AutoIntegral_lag = zeros(NCells, NTimeLag, NAutocorrFn);
lag_list = [0:(NTimeLag-1)]' * SamplingInterval;

for cell_k = 1:NCells
    % 0 lag: store EKE
    t_end = 1;
    AutoIntegral_lag(cell_k, t_end, :) = AutocorrFn(cell_k, t_end, :);
    
    % +ve lag: store diffusivity
    for t_end = 2:NTimeLag
        lag_clipped = lag_list(1:t_end);
        AutocorrFn_clipped = squeeze(AutocorrFn(cell_k, 1:t_end, :));
        
        AutoIntegral_lag(cell_k, t_end, :) = trapz(lag_clipped, AutocorrFn_clipped);
    end
end

end