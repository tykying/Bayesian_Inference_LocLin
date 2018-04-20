function [TrajJumps_DA] = obtain_TrajJumps_Ito_DataAnalysis(TrajJumps, Mesh)
    % Compute Data Statistics(Moments)
    DataMoments = zeros(1, 4);
    TrajJumps_Stat = struct();
    
    TrajJumps_MomentGlobal = struct();
    % Initialise TrajJumps_StatGlobal
    cell_I = 1;
    fields = fieldnames(TrajJumps(cell_I));
    for i = 1:numel(fields)
        TrajJumps_MomentGlobal.(fields{i}) =  zeros(length(TrajJumps), 4);
    end
    
    % Fill in the Data
    for cell_I = 1:length(TrajJumps)
        fields = fieldnames(TrajJumps(cell_I));
        for i = 1:numel(fields)
          Data = TrajJumps(cell_I).(fields{i});
          
          DataMoments(1) = mean(Data);
          DataMoments(2) = var(Data);
          DataMoments(3) = skewness(Data);
          DataMoments(4) = kurtosis(Data);
          
          TrajJumps_Stat(cell_I).(fields{i}) =  DataMoments;  % Not advise to use
          TrajJumps_MomentGlobal.(fields{i})(cell_I, :) = DataMoments;
        end
    end
    disp('obtain_TrajJumps_DataAnalysis: Finished Computing JumpsCFL.')

    % Compute CFL number of Jumps
    JumpsCFL = zeros(length(TrajJumps), 2);
    for cell_I = 1:length(TrajJumps)
        diffx = TrajJumps(cell_I).diffx;
        diffy = TrajJumps(cell_I).diffy;
        
        h = TrajJumps(cell_I).h;
        JumpsCFL(cell_I, 1) = mean(abs(diffx))/Mesh(cell_I).cell_dx(1);
        JumpsCFL(cell_I, 2) = mean(abs(diffy))/Mesh(cell_I).cell_dx(2);
    end
    
    disp('obtain_TrajJumps_DataAnalysis: Finished Computing JumpsCFL.')
    TrajJumps_DA = struct('JumpsCFL', JumpsCFL, 'TrajJumps_MomentGlobal', TrajJumps_MomentGlobal);
    
    disp('obtain_TrajJumps_DataAnalysis: Finished All Jumps Data Analysis.')
end