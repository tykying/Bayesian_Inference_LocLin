function [Data_ijk_FILTERED] = mean_filter_Dataijk(Data_ijk)
    Data_ijk_FILTERED = Data_ijk;
    
    for k = 1:size(Data_ijk, 3)
        for j = 1:size(Data_ijk, 2)
            for i = 1:size(Data_ijk, 1)
                i_bgn = max(1, i-1); i_end = min(i+1, size(Data_ijk, 1));
                j_bgn = max(1, j-1); j_end = min(j+1, size(Data_ijk, 2));
                
                Data_ijk_BLOCK = Data_ijk(i_bgn:i_end, j_bgn:j_end, k);
                Data_ijk_FILTERED(i, j, k) = mean(Data_ijk_BLOCK(:));
            end
        end
    end
    
end
    