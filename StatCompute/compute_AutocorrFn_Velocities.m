function [AutocorrFn] = compute_AutocorrFn_Velocities(Autocorr_Sum, AvailSample)
    [NCells, NTimeLag, NAutocorrFn] = size(Autocorr_Sum);
    
    AutocorrFn = zeros(size(Autocorr_Sum));
    for cell_k = 1:NCells
        AvailSample_cell = AvailSample(cell_k, :);
        AutocorrFn(cell_k, :, :) = squeeze(Autocorr_Sum(cell_k, :, :))./AvailSample(cell_k, :)';
        
        % AutocorrFn(cell_k, :, :) = AutocorrFn(cell_k, :, :)./AutocorrFn(cell_k, 1, :);  %% Normalise w.r.t. lag=0 day
    end
end