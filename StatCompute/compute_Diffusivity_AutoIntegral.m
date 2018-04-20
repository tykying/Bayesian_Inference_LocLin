function [K, K_Conv_sdr] = compute_Diffusivity_AutoIntegral(AutoIntegral, AvailSample, min_samplesize)

[NCells, NTimeLag, NAutocorrFn] = size(AutoIntegral);

K = zeros(NCells, NAutocorrFn);
K_Conv_sdr = zeros(NCells, NAutocorrFn);

for cell_k = 1:NCells
    tau_filter = (AvailSample(cell_k, :) > min_samplesize);
    tau_filter_len = sum(tau_filter);
    
    tau_clipped = [floor(tau_filter_len*0.9):tau_filter_len];
    AutoIntegral_clipped = squeeze(AutoIntegral(cell_k, tau_clipped, :));
    
    % Check if convergence is reached
    K(cell_k, :) = mean(AutoIntegral_clipped);
    K_Conv_sdr(cell_k, :) = sqrt(var(AutoIntegral_clipped))./abs(K(cell_k, :));
end

end
