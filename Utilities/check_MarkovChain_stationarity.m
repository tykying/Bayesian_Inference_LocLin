% % Good 
% MC_Data = theta_store(79, 1, :);
% 
% % Bad
% MC_Data = theta_store(235, 1, :);


function MCStat_TF = check_MarkovChain_stationarity(MC_Data)
%% MC Stationary
Npartitions = 4;
assert(Npartitions >= 3)

Mean_section = zeros(1,Npartitions);
Vari_section = zeros(1,Npartitions);

NSamples = length(MC_Data);

for partition = 1:Npartitions
    partition_bgn = round(NSamples*(partition-1)/Npartitions)+1;
    partition_end = round(NSamples*partition/Npartitions);
    %fprintf('Partition: %d, in samples (%d, %d). \n', partition, partition_bgn, partition_end);
    
    MC_Data_partition = MC_Data(partition_bgn:partition_end);
    
    Mean_section(partition) = mean(MC_Data_partition);
    Vari_section(partition) = var(MC_Data_partition);
end

%% Compare the final section with the previous two sectionss
% Assume the Final Three are all correct distribution -
% Test against each other whether they have the same mean
Nref = 3;
MCStat_MatrixTF = ones(Nref, Nref);
z_c = 1.64;  % 90% CI

ref_ind = 0;
for ref = (Npartitions-(Nref-1)):Npartitions
    ref_ind = ref_ind + 1;
    Mean_ref = Mean_section(ref);
    Vari_ref = Vari_section(ref);
    
    sec_ind = 0;
    for sec = (Npartitions-(Nref-1)):Npartitions
        sec_ind = sec_ind + 1;
        z_test = abs((Mean_section(sec) - Mean_ref)/sqrt(Vari_ref));
        
        if z_test > z_c
            MCStat_MatrixTF(ref_ind, sec_ind) = 0;
        end
    end
end

MCStat_TF = all(MCStat_MatrixTF(:)==1);

end