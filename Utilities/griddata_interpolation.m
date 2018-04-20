function [refined_data] = griddata_interpolation(data, refinement_ratio)
    % refinement_ratio: row vector of two elements; Ratio to refine
    input_data_size = size(data);
    
    refined_data = zeros(input_data_size.*refinement_ratio);
    
    for i_coarse = 1:input_data_size(1)
        for j_coarse = 1:input_data_size(2)
            
            for local_ri = (refinement_ratio(1)-1):-1:0
            	for local_rj = (refinement_ratio(2)-1):-1:0
                    i_fine = refinement_ratio(1)*i_coarse-local_ri;
                    j_fine = refinement_ratio(2)*j_coarse-local_rj;
                    refined_data(i_fine, j_fine) = data(i_coarse, j_coarse);
                end
            end
        end
    end
    
    
end
