function [theta_MAP, logPost_MAP, MAP_SampleInd] = locate_theta_MAP(theta_store, logPost_ts)

Ncells = size(theta_store, 1);
NVars = size(theta_store, 2);

logPost_MAX = max(logPost_ts, [], 2);
theta_MAP = zeros(Ncells, NVars);
MAP_SampleInd = zeros(Ncells, 1);
logPost_MAP = zeros(Ncells, 1);

for cell_k = 1:Ncells
    MAP_SampleInd(cell_k) = find( logPost_ts(cell_k, :) == logPost_MAX(cell_k), 1);
    
    logPost_MAP(cell_k) = logPost_MAX(cell_k);
    theta_MAP(cell_k, :) = squeeze( theta_store(cell_k,:,MAP_SampleInd(cell_k)) );
end


end