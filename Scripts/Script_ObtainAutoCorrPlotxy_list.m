% Given cell_k_list, min_samplesize
plotx_list = cell(size(cell_k_list));
ploty_list = cell(size(cell_k_list));

for list_iterator = 1:length(cell_k_list(:))
    cell_k = cell_k_list(list_iterator);
    
    tau_filter = (AvailSample(cell_k, :) > min_samplesize);
    tau_len = sum(tau_filter);
    
    time_lag = [0:(tau_len-1)]*SamplingInterval/time_scale;  % starting from lag=0
    
    plotx_data = time_lag';
    ploty_data = squeeze(AutocorrFn(cell_k, tau_filter, CPN_vis));
    
    plotx_list{list_iterator} = plotx_data;
    ploty_list{list_iterator} = ploty_data;
end