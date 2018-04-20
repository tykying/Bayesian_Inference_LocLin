function [TransMat] = obtain_TransMat(TrajData, TrajPosition, Mesh)

% Optimised; But not directly parallelisable (due to data traffic)
[TrajJumps_iBased, NJumps_cell] = sort_TrajJumps_iBased(TrajData, TrajPosition, Mesh);

% Formulate Transition Matrix: P(i, j) = P(j|i)
TransMat = zeros(length(NJumps_cell), length(NJumps_cell));
for cell_I = 1:length(TrajJumps_iBased)
    TransDen_rowi = zeros(1, length(NJumps_cell));
    for jump_celli = 1:NJumps_cell(cell_I)
        cellj_F = TrajJumps_iBased(cell_I).GridInd_f(jump_celli);
        if (cellj_F >= 1) && (cellj_F <= length(Mesh))
            TransDen_rowi(cellj_F) = TransDen_rowi(cellj_F) + 1;
        end
    end
    
    TransMat(cell_I, :) = TransDen_rowi;
end

end