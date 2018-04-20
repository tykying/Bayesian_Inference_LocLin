% Store the initial and destinated alpha

% Reminder:
% Group in destinated cell (fBased) seems more efficient to memory use: only need to
% load the different realisations of LOCAL likelihood function to perform interpolations
function [TrajJumps, NJumps_cell] = sort_TrajJumps_iBased(TrajData, TrajPosition, Mesh)
    NCells = length(Mesh);
    NJumps_cell = zeros(1, NCells);
    TrajJumps(NCells) = struct('alphax_i', [], 'alphay_i', [], ...
                            'alphax_f', [], 'alphay_f', [], ...
                            'h', [], 'GridInd_f', []);
%                           'GridInd_i', []);,  ...

    % For iBased: no need to store GridInd_i!
    
    % Assumed uniform sampling interval
    ts_list = TrajData.ts_list;
    
    GridInd = TrajPosition.GridInd;
    alphax = TrajPosition.alphax;
    alphay = TrajPosition.alphay;
    
    [Nts, nparticles] = size(GridInd);
        
    % Assume the memory is big enough to store all _i, _f
    GridInd_i = GridInd(1:(Nts-1), :);
    GridInd_f = GridInd(2:Nts, :);
    alphax_i = alphax(1:(Nts-1), :);
    alphax_f = alphax(2:Nts, :);
    alphay_i = alphay(1:(Nts-1), :);
    alphay_f = alphay(2:Nts, :);
    
    % Make the h_list a full matrix rather than an array
    h_list = diff(ts_list).*ones(Nts-1, nparticles);
    assert(all(size(h_list) == size(GridInd_i)));
    
    for cell_k = 1:NCells  % Parallelise: each node takes care of some cells
        if (mod(log2(cell_k), 1) == 0)
            disp(['Sorting Jumps in cell', num2str(cell_k)]);
        end
        
        % Initialise TrajJumps
        TrajJumps_cellk = TrajJumps(cell_k);
                
        GridInd_filter_cellk = (GridInd_i == cell_k);  % i_Based
        
        TrajJumps_cellk.alphax_i = alphax_i(GridInd_filter_cellk);
        TrajJumps_cellk.alphay_i = alphay_i(GridInd_filter_cellk);
        
        TrajJumps_cellk.GridInd_f = GridInd_f(GridInd_filter_cellk);
        TrajJumps_cellk.alphax_f = alphax_f(GridInd_filter_cellk);
        TrajJumps_cellk.alphay_f = alphay_f(GridInd_filter_cellk);

        TrajJumps_cellk.h = h_list(GridInd_filter_cellk);
        
        TrajJumps(cell_k) = TrajJumps_cellk;
        NJumps_cell(cell_k) = length(TrajJumps_cellk.GridInd_f);
    end
    
end