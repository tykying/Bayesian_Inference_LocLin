function SampleMomentsStruct = evaluate_SampleMomentsStruct(SampleData)

SampleDataFields = fieldnames(SampleData(1));
N_Fields = numel(SampleDataFields);

N_Samples = numel(SampleData);
Ncell_SpDis = size(SampleData(1).(SampleDataFields{1}), 1);


% Initialise SampleMomentsStruct
SampleMomentsStruct = struct;
for i = 1:numel(SampleDataFields)   
    SampleMomentsStruct.(SampleDataFields{i}).Mean = zeros(size(SampleData(1).(SampleDataFields{i})));
    SampleMomentsStruct.(SampleDataFields{i}).Vari = zeros(size(SampleData(1).(SampleDataFields{i})));
    SampleMomentsStruct.(SampleDataFields{i}).Skew = zeros(size(SampleData(1).(SampleDataFields{i})));
    SampleMomentsStruct.(SampleDataFields{i}).Kurt = zeros(size(SampleData(1).(SampleDataFields{i})));    
end

% Compute SampleMomentsStruct
for cell = 1 : Ncell_SpDis
    for i = 1 : N_Fields
        % Obtain the correct data structure for the Data
        [RawDataNCell, RawDataDim] = size(SampleData(1).(SampleDataFields{i}));  % First dimension = number of cells
        
        CellRawData = zeros(N_Samples, RawDataDim);
        for Sample = 1:N_Samples
            CellRawData(Sample, :) = SampleData(Sample).(SampleDataFields{i})(cell, :);
        end
        
        SampleMomentsStruct.(SampleDataFields{i}).Mean(cell,:) = mean(CellRawData);
        SampleMomentsStruct.(SampleDataFields{i}).Vari(cell,:) = var(CellRawData);
        SampleMomentsStruct.(SampleDataFields{i}).Skew(cell,:) = skewness(CellRawData);
        SampleMomentsStruct.(SampleDataFields{i}).Kurt(cell,:) = kurtosis(CellRawData);        
    end
end

end
