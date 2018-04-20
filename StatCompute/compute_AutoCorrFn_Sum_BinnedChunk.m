function [Autocorr_Sum_BinnedChunk, AvailSample_BinnedChunk] = compute_AutoCorrFn_Sum_BinnedChunk(TrajData, TrajPosition, VelocData, TransitStat_Part, Mesh_Struct)
%% Explanation:
% BinnedChunk: consider segments of trajectories fully in the bin

NCells = length(Mesh_Struct.Mesh);

x = TrajData.x;
y = TrajData.y;
ts_list = TrajData.ts_list;

GridInd = TrajPosition.GridInd;

u = VelocData.u;
v = VelocData.v;
ds = sqrt(u.^2 + v.^2);

[Nts, nparticles] = size(x);

% 5 Components: XX, YY, XY, YX, dsds
NCpn = 5;
BinnedChunk_MaxLength = 365;
Autocorr_Sum_BinnedChunk = zeros(NCells, BinnedChunk_MaxLength, NCpn);
AvailSample_BinnedChunk = zeros(NCells, BinnedChunk_MaxLength);

Autocorr_BinnedChunk_cellk = zeros(1, BinnedChunk_MaxLength, NCpn);

for part = 1:nparticles
    ArrivalTime = TransitStat_Part(part).ArrivalTime;
    ResideIntv = TransitStat_Part(part).ResideIntv;

    % Number of available segments = Number of hitting
    AvailChunk = length(ArrivalTime);
    
    % TempChunk: consider chunks of trajectories fully in the bin
    for Chunk_ind = 1:AvailChunk
        % Ignore the statistics if the residence time is beyond MaxLength
        ChunkBinResideIntv = min(ResideIntv(Chunk_ind), BinnedChunk_MaxLength);
        
        % The beginning and ending i_t of the Chunk
        Chunk_i_ti = ArrivalTime(Chunk_ind);
        Chunk_i_tf = ArrivalTime(Chunk_ind) + (ChunkBinResideIntv-1);  % In the time-series of position, not diff(position)!
        
        % During this interval do not wish GridInd to vary too much...
        Chunk_X_ind = Chunk_i_ti:Chunk_i_tf;
        
        Chunk_GridInd = GridInd(Chunk_i_ti, part);
%         % For debugging; Should be removed for efficiency
%         assert(length(Chunk_X_ind) == ((ChunkBinResideIntv-1) + 1));
%         assert(Chunk_GridInd == unique(GridInd(Chunk_X_ind, part)));
%         if (ChunkBinResideIntv ~= BinnedChunk_MaxLength)
%             if (Chunk_i_tf ~= Nts)
%                 assert(Chunk_GridInd ~= GridInd(Chunk_i_tf+1, part));
%             end
%         end
        
        % Select the jumps within the chunk interval
        Chunk_JumpsLen = ChunkBinResideIntv-1;
        Chunk_U_ind = Chunk_i_ti:(Chunk_i_tf-1);  % -1: to change from the view point of position to diff(position)
        
        Chunk_u = u(Chunk_U_ind, part);
        Chunk_v = v(Chunk_U_ind, part);
        Chunk_ds = ds(Chunk_U_ind, part);
        
        %assert(length(Chunk_U_ind) == (ChunkBinResideIntv-1));
        
        % Cross-Correlation between time-series 
        % If ts1==ts2: Autocorrelation of a single time-series
        % If ts1~=ts2: Mix of auto and cross-correlation; also XY~=YX
        LagInd = (0:(Chunk_JumpsLen-1))+1;
        CrossSum_NSample = [Chunk_JumpsLen:-1:1];
        
        AvailSample_BinnedChunk(Chunk_GridInd, LagInd) = AvailSample_BinnedChunk(Chunk_GridInd, LagInd) + CrossSum_NSample;
        
        % Eddy Correlation
        Autocorr_BinnedChunk_cellk(1, LagInd, 1) = compute_LagCrossCorrSum(Chunk_u, Chunk_u);
        Autocorr_BinnedChunk_cellk(1, LagInd, 2) = compute_LagCrossCorrSum(Chunk_v, Chunk_v);
        Autocorr_BinnedChunk_cellk(1, LagInd, 3) = compute_LagCrossCorrSum(Chunk_u, Chunk_v);
        Autocorr_BinnedChunk_cellk(1, LagInd, 4) = compute_LagCrossCorrSum(Chunk_v, Chunk_u);
        Autocorr_BinnedChunk_cellk(1, LagInd, 5) = compute_LagCrossCorrSum(Chunk_ds, Chunk_ds);

        Autocorr_Sum_BinnedChunk(Chunk_GridInd, LagInd, :) = Autocorr_Sum_BinnedChunk(Chunk_GridInd, LagInd, :) + Autocorr_BinnedChunk_cellk(1, LagInd, :);
        %Autocorr_BinnedChunk_cellk(:) = 0;
    end
    
    % Report Progress
    if mod(log2(part), 2) == 0
        disp(['Finished particle ', num2str(part)]);
    end
end

end