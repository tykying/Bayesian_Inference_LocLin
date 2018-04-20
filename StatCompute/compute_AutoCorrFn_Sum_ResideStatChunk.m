function [Autocorr_Sum_TemporChunk, AvailSample_TemporChunk] = compute_AutoCorrFn_Sum_ResideStatChunk(TrajData, TrajPosition, VelocData, TransitStat_Part, Mesh_Struct, ResideIntv_Stat)
%% Explanation:
% TemporChunk: consider segments of trajectories fully in the bin

NCells = length(Mesh_Struct.Mesh);
Nx_cell = Mesh_Struct.RectMesh_Param.Nx_cell;

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
ChunkBinResideIntv_min = 2;
TemporChunk_MaxLength = 240;
Autocorr_Sum_TemporChunk = zeros(NCells, TemporChunk_MaxLength, NCpn);
AvailSample_TemporChunk = zeros(NCells, TemporChunk_MaxLength);

Autocorr_TemporChunk_cellk = zeros(1, TemporChunk_MaxLength, NCpn);

for part = 1:nparticles
    ArrivalTime = TransitStat_Part(part).ArrivalTime;
    ResideIntv = TransitStat_Part(part).ResideIntv;

    % Only consider trajectory segments that stay in its initial cell for
    % at least some time-steps
    IndepChunk_filter = (ResideIntv > ChunkBinResideIntv_min);
    
    % Number of available segments = Number of hitting
    ArrivalTime = ArrivalTime(IndepChunk_filter);
    %ResideIntv = ResideIntv(IndepChunk_filter);
    
    AvailChunk = length(ArrivalTime);
    
    % TempChunk: consider a fixed length of chunks of trajectories fully in the bin
    for Chunk_ind = 1:AvailChunk      
        % The beginning and ending i_t of the Chunk
        Chunk_i_ti = ArrivalTime(Chunk_ind);
        Chunk_GridInd = GridInd(Chunk_i_ti, part);

        %% Adaptation here:
        % Select length of chunk by statistics of ResidenceIntv
        ResideChunk_MaxLength = 2*ResideIntv_Stat(Chunk_GridInd, 3)*(Nx_cell/8);  % Length = 2*Median*(Nx_Cell/8)
        Chunk_i_tf = min([Chunk_i_ti+ResideChunk_MaxLength, Chunk_i_ti+TemporChunk_MaxLength, Nts]);  % In the time-series of position, not diff(position)!
        % Avoid the chunk exceeds the time limit
        
        % During this interval do not wish GridInd to vary too much...
        Chunk_X_ind = Chunk_i_ti:Chunk_i_tf;
        ChunkBinResideIntv = length(Chunk_X_ind);        
                
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
        
        AvailSample_TemporChunk(Chunk_GridInd, LagInd) = AvailSample_TemporChunk(Chunk_GridInd, LagInd) + CrossSum_NSample;
        
        % Eddy Correlation
        Autocorr_TemporChunk_cellk(1, LagInd, 1) = compute_LagCrossCorrSum(Chunk_u, Chunk_u);
        Autocorr_TemporChunk_cellk(1, LagInd, 2) = compute_LagCrossCorrSum(Chunk_v, Chunk_v);
        Autocorr_TemporChunk_cellk(1, LagInd, 3) = compute_LagCrossCorrSum(Chunk_u, Chunk_v);
        Autocorr_TemporChunk_cellk(1, LagInd, 4) = compute_LagCrossCorrSum(Chunk_v, Chunk_u);
        Autocorr_TemporChunk_cellk(1, LagInd, 5) = compute_LagCrossCorrSum(Chunk_ds, Chunk_ds);

        Autocorr_Sum_TemporChunk(Chunk_GridInd, LagInd, :) = Autocorr_Sum_TemporChunk(Chunk_GridInd, LagInd, :) + Autocorr_TemporChunk_cellk(1, LagInd, :);
        %Autocorr_TemporChunk_cellk(:) = 0;
    end
    
    % Report Progress
    if mod(log2(part), 2) == 0
        disp(['Finished particle ', num2str(part)]);
    end
end

end