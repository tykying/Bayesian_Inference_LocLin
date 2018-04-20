% Store the initial and destinated alpha

% Reminder:
% Group in destinated cell (fBased) seems more efficient to memory use: only need to
% load the different realisations of LOCAL likelihood function to perform interpolations
function [TrajJumps, TrajJumps_Length] = sort_TrajJumps_fBased(TrajData, TrajPosition, Mesh)
    NCells = length(Mesh);
    TrajJumps_Length = zeros(1, length(Mesh));

    
    % Assumed uniform sampling interval
    ts_list = TrajData.ts_list;
    
    GridInd = TrajPosition.GridInd;
    alphax = TrajPosition.alphax;
    alphay = TrajPosition.alphay;
    
    [Nts, nparticles] = size(GridInd);
    
    % Loop over all data: 
    % Count size of array needed
    for cell_k = 1:NCells  % Parallelise: each node takes care of some cells
        for t_ind = 2:length(ts_list)  % fBased; iBased -> t_ind = 1:(length(ts_list)-1)
            GridInd_t_ind = GridInd(t_ind, :);
            GridInd_filter = (GridInd_t_ind == cell_k);
            
            TrajJumps_Length(cell_k) = TrajJumps_Length(cell_k) + sum(GridInd_filter);
        end
    end
    
    for cell_k = 1:NCells  % Parallelise: each node takes care of some cells
        cell_Length = TrajJumps_Length(cell_k);
        
        % Initialise TrajJumps
        TrajJumps(cell_k) = struct('alphax_i', zeros(1, cell_Length),  ...
                            'alphay_i', zeros(1, cell_Length),  ...
                            'alphax_f', zeros(1, cell_Length),  ...
                            'alphay_f', zeros(1, cell_Length),  ...
                            'GridInd_i', zeros(1, cell_Length),  ...
                            'GridInd_f', zeros(1, cell_Length),  ...
                            'h', zeros(1, cell_Length));
    end
    
    for cell_k = 1:NCells  % Parallelise: each node takes care of some cells
        Length_counter = 0;
        
        for t_ind = 2:length(ts_list)  % fBased
            for part = 1:nparticles
                if (GridInd(t_ind, part)==cell_k)
                    Length_counter = Length_counter + 1;

                    % fBased; iBased -> take t_ind, t_ind+1
                    TrajJumps(cell_k).GridInd_i(Length_counter) = GridInd(t_ind-1, part);
                    TrajJumps(cell_k).alphax_i(Length_counter) = alphax(t_ind-1, part);
                    TrajJumps(cell_k).alphay_i(Length_counter) = alphay(t_ind-1, part);
                    
                    TrajJumps(cell_k).GridInd_f(Length_counter) = GridInd(t_ind, part);
                    TrajJumps(cell_k).alphax_f(Length_counter) = alphax(t_ind, part);
                    TrajJumps(cell_k).alphay_f(Length_counter) = alphay(t_ind, part);
                    
                    TrajJumps(cell_k).h(Length_counter) = ts_list(t_ind)-ts_list(t_ind-1);
                end
            end
        end
        
        assert(Length_counter == TrajJumps_Length(cell_k));
    end
    
end