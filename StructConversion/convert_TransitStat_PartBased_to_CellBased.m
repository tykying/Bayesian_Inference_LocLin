function [TransitStat_Cell, EffectiveNTraj_Cell] = convert_TransitStat_PartBased_to_CellBased(TransitStat_Part, Mesh_Struct)
%% Term Explained 
% ArrivalTime: time at which the particle transits to a new cell (in terms of i_t)
% ResideIntv: time interval for a particle before the next transition (in terms of i_t)

% Time line:
% time:     |-t1-||-t2-||-t3-||-t4-|
% absPos:   |-x1-||-x2-||-x3-||-x4-|
% diffPos:  X----u1----u2----u3----X

% CellInd:  |-C1-||-C1-||-C2-||-C3-|
% -> ResideIntv in C1 = 2; ArrivalTime at C1 = 1;
% -> ResideIntv in C2 = 1; ArrivalTime at C2 = 3;

% Remark:
% Initial position is COUNTED as the first arrival
%% Computation
NCells = length(Mesh_Struct.Mesh);
nparticles = length(TransitStat_Part);

%   Initialise the struct for output
EffectiveNTraj_Cell = zeros(NCells, 3);

%   Predetermine the sizes of arrays
array_size_I = zeros(NCells, 1); % For 'ResideIntv_I', 'GridInd_To_I' and 'Part_ID_I'
array_size_Itm = zeros(NCells, 1);
array_size_F = zeros(NCells, 1);

for part = 1: nparticles
    Part_GridInd_List = TransitStat_Part(part).GridInd_List;

    % Initial cell
    SampleInd = 1;
    
    cell_k = Part_GridInd_List(SampleInd);
    array_size_I(cell_k) = array_size_I(cell_k)+1;
    
    % Intermediate cells; Exclude the final cell
    for SampleInd = 2:(length(Part_GridInd_List)-1)
        cell_k = Part_GridInd_List(SampleInd);
        array_size_Itm(cell_k) = array_size_Itm(cell_k)+1;
    end
    
    % Final cell
    SampleInd = length(Part_GridInd_List);
    
    cell_k = Part_GridInd_List(SampleInd);
    array_size_F(cell_k) = array_size_F(cell_k)+1;
end

%   Initialise the struct for output
TransitStat_Cell = struct('ResideIntv', [], ...
    'ResideIntv_I', [], 'ResideIntv_F', [], ...
    'GridInd_To_I', [], 'GridInd_From_F', [], ...
    'GridInd_From', [], 'GridInd_To', [], ...
    'ArrivalTime', [], 'ArrivalTime_F', [], ...
    'Part_ID', [], ...
    'Part_ID_I', [], 'Part_ID_F', []);

for cell_k = 1:NCells
    size_I = [array_size_I(cell_k), 1];
    size_Itm = [array_size_Itm(cell_k), 1];
    size_F = [array_size_F(cell_k), 1];

    TransitStat_Cell(cell_k) = struct('ResideIntv', zeros(size_Itm), ...
    'ResideIntv_I', zeros(size_I), 'ResideIntv_F', zeros(size_F), ...
    'GridInd_To_I', zeros(size_I), 'GridInd_From_F', zeros(size_F), ...
    'GridInd_From', zeros(size_Itm), 'GridInd_To', zeros(size_Itm), ...
    'ArrivalTime', zeros(size_Itm), 'ArrivalTime_F', zeros(size_F), ...
    'Part_ID', zeros(size_Itm), ...
    'Part_ID_I', zeros(size_I), 'Part_ID_F', zeros(size_F));
end

%   Fill up the arrays
i_I = zeros(NCells, 1);
i_Itm = zeros(NCells, 1);
i_F = zeros(NCells, 1);

for part = 1: nparticles
    % Load PartBased data
    Part_ResideIntv = TransitStat_Part(part).ResideIntv;
    Part_ArrivalTime = TransitStat_Part(part).ArrivalTime;
    Part_GridInd_List = TransitStat_Part(part).GridInd_List;
    
    % Initial cell
    SampleInd = 1;
    
    cell_k = Part_GridInd_List(SampleInd);
    i_I(cell_k) = i_I(cell_k) + 1;   % Increment counter in each cell
    ind = i_I(cell_k);
    
    Part_GridInd_List_len = length(Part_GridInd_List);
    
    TransitStat_Cell(cell_k).ResideIntv_I(ind) = Part_ResideIntv(SampleInd);
    
    % Special treatment for particles which remain in one cell
    if (Part_GridInd_List_len > 1)
        TransitStat_Cell(cell_k).GridInd_To_I(ind) = Part_GridInd_List(SampleInd+1);
    else
        TransitStat_Cell(cell_k).GridInd_To_I(ind) = 0;
    end

    TransitStat_Cell(cell_k).Part_ID_I(ind) = part;

    % Intermediate cells; Exclude the final cell
    for SampleInd = 2:(Part_GridInd_List_len-1)
        cell_k = Part_GridInd_List(SampleInd);
        i_Itm(cell_k) = i_Itm(cell_k) + 1;
        ind = i_Itm(cell_k);
        
        TransitStat_Cell(cell_k).ResideIntv(ind) = Part_ResideIntv(SampleInd);
        TransitStat_Cell(cell_k).ArrivalTime(ind) = Part_ArrivalTime(SampleInd);

        TransitStat_Cell(cell_k).GridInd_From(ind) = Part_GridInd_List(SampleInd-1);
        TransitStat_Cell(cell_k).GridInd_To(ind) = Part_GridInd_List(SampleInd+1);
        
        TransitStat_Cell(cell_k).Part_ID(ind) = part;
    end
    
    % Final cell
    SampleInd = Part_GridInd_List_len;
    
    cell_k = Part_GridInd_List(SampleInd);
    i_F(cell_k) = i_F(cell_k) + 1;
    ind = i_F(cell_k);

    TransitStat_Cell(cell_k).ResideIntv_F(ind) = Part_ResideIntv(SampleInd);
    TransitStat_Cell(cell_k).ArrivalTime_F(ind) = Part_ArrivalTime(SampleInd);
    
    % Special treatment for particles which remain in one cell
    if (Part_GridInd_List_len > 1)
        TransitStat_Cell(cell_k).GridInd_From_F(ind) = Part_GridInd_List(SampleInd-1);
    else
        TransitStat_Cell(cell_k).GridInd_From_F(ind) = 0;
    end

        
    TransitStat_Cell(cell_k).Part_ID_F(ind) = part;
end

% Debugging
assert(all(i_I == array_size_I));
assert(all(i_Itm == array_size_Itm));
assert(all(i_F == array_size_F));

% Effective Trajectories in each cell
EffectiveNTraj_Cell = [i_I, i_Itm, i_F];
end
