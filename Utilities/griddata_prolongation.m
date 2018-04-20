function [coarsened_data] = griddata_prolongation(data, coarsening_ratio)
    input_data_size = size(data);
    
    assert(all(coarsening_ratio >= 1))
    
    coarsened_data = zeros(input_data_size./coarsening_ratio);
    
    for i_coarse = 1:size(coarsened_data, 1)
        for j_coarse = 1:size(coarsened_data, 2)
            
            local_ri_range = [coarsening_ratio(1)*(i_coarse-1)+1, coarsening_ratio(1)*i_coarse];
            local_rj_range = [coarsening_ratio(2)*(j_coarse-1)+1, coarsening_ratio(2)*j_coarse];
            
            coarsened_data(i_coarse, j_coarse) = mean(mean(data(local_ri_range, local_rj_range)));
        end
    end
    
    
end
