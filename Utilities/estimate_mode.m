function Peak_est = estimate_mode(SamplesStored)
MaxNPeaks = 7;
nbins_List = 8:1:64;
edge_c_Peaks = zeros(length(nbins_List), MaxNPeaks);

% Records of peaks under different resolutions
for nbins_ind = 1:length(nbins_List)
    nbins = nbins_List(nbins_ind);
    
    [N,edges] = histcounts(SamplesStored, nbins, 'Normalization', 'pdf');
    edges_c = 0.5*( edges(1:end-1) + edges(2:end) );
    % pks =findpeaks(N)
    
    [N_val N_ind] = sort(N, 'descend');
    edge_c_Peaks(nbins_ind, :) = edges_c(N_ind(1:MaxNPeaks));
end

% Coarse histogram for the locations of the peaks
[N_Peaks,edges_Peak] = histcounts(edge_c_Peaks(:), 16);
max_ind =  find(N_Peaks == max(N_Peaks));

% Select a region in which the mean will be computed
FilterRegion = [edges_Peak(max_ind), edges_Peak(max_ind+1)];
FilterInd = ( (SamplesStored > FilterRegion(1)) ...
    & (SamplesStored < FilterRegion(2)) );
SamplesFiltered = SamplesStored(FilterInd);

Peak_est = mean(SamplesFiltered);
end