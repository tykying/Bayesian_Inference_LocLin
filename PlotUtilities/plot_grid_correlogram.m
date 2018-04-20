% x, y, ts_list
ts_list = linspace(0, 10*pi,1000)';
x = [0.95*sin(3*ts_list), 0.90*sin(8*(ts_list+0.5))];
y = [0.90*cos(1*ts_list), 0.95*cos(12*(ts_list+0.5))];
%plot(x,y)

L = 1;

% Use saved statistics
Use_Saved_Stat = 0;

% TS Tracking Inference
TrajType_Save = 'ResidTime_QGM2';
%DataTypeString = '_FrequStat';
%nparticles = 2000;
Nx_cell_options = 32;
DGscheme = 'U0K0';

%load('/home/s1046972/opt/qgm2_particle/PART_TRAJ_nu300/traj_Npart10000_tIntv5yr_SampletIntv24hr.mat')
%load('/home/s1046972/opt/qgm2_particle/PART_TRAJ/traj_Npart2000_tIntv9yr_SampletIntv24hr.mat')
load('/home/s1046972/opt/qgm2_particle/PART_TRAJ/traj_Npart20000_tIntv8yr_SampletIntv24hr.mat')

% Read x,y,ts_list,L,ngrid
%%% Convert x, y from centimeters to meters
x = x/(100); y = y/(100);  L = L/(100);

%h = (ts_list(2)-ts_list(1))/(3600);
h = (ts_list(2)-ts_list(1));


% Particle Location
SamplingInterval = 1*3600*24;
SamplingParticle = 1;

SpDis_ranges_min = [0, 0]; SpDis_ranges_max = [L, L];
Nx_cell = Nx_cell_options; Ny_cell = Nx_cell_options;

[Mesh, RectMesh_Param] = setup_RectMeshStruct(Nx_cell, Ny_cell, SpDis_ranges_min, SpDis_ranges_max);
Mesh_Struct = struct('Mesh', Mesh, 'RectMesh_Param', RectMesh_Param);

grid_x = RectMesh_Param.grid_x;
grid_y = RectMesh_Param.grid_y;

cellij_to_k = RectMesh_Param.cellij_to_k;
cellk_to_ij = RectMesh_Param.cellk_to_ij;


%load('./Trajectories_Generated/taylor_green_with_shear_dudy0.1_const_Pe10_trajectories_saved.mat')
%stop

% Some interesting Trajectories
% SamplingInterval = 1*3600*24; SamplingParticle = 250: curling in subpolar gyre

TrajData_SubSampled = SubSample_Traj(x, y, ts_list, SamplingInterval, SamplingParticle);

%plot(TrajData_SubSampled.x_SubSampled, TrajData_SubSampled.y_SubSampled)

SamplingInterval_ind = SamplingInterval/h;

space_scale = 1000;
time_scale = 3600;
SamplingInterval_vis = SamplingInterval/time_scale;


nparticles = size(TrajData_SubSampled.x_SubSampled, 2);
Nts = length(TrajData_SubSampled.ts_list_SubSampled);


DataTypeString = '_TransStat';
filename = ['/home/s1046972/data/', TrajType_Save, DataTypeString, '_npart', num2str(nparticles), '_SamplingIntv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell) '.mat'];

if Use_Saved_Stat == 0
  
    GridInd = zeros(Nts, nparticles);
    
    alphax = zeros(Nts, nparticles);
    alphay = zeros(Nts, nparticles);
        
    % Locate the particle in each bin
    for i_t = 1: Nts        
        for i_part = 1: nparticles
            x_part = x(i_t, i_part);
            y_part = y(i_t, i_part);
            
            [CellIndi_part, alphax_part] = search_sorted_nonuniform(x_part, grid_x);   %%%%% HARDCODED: 2D Field here
            [CellIndj_part, alphay_part] = search_sorted_nonuniform(y_part, grid_y);   %%%%% HARDCODED: 2D Field here
            
            CellIndk_part = cellij_to_k(CellIndi_part, CellIndj_part);
            
            GridInd(i_t, i_part) = CellIndk_part;
            
            alphax(i_t, i_part) = alphax_part;
            alphay(i_t, i_part) = alphay_part;
        end
    end
    
    Analysed_Traj = struct('x', TrajData_SubSampled.x_SubSampled, ...
                       'y', TrajData_SubSampled.y_SubSampled, ...
                       'ts_list', TrajData_SubSampled.ts_list_SubSampled, ...
                       'GridInd', GridInd, ...
                       'alphax', alphax, ...
                       'alphay', alphay )
    
    % Statistics of Arrival time and Residence time
    diff_GridInd = diff(GridInd, 1);
    
    Transit_Part = struct('Arrival_time', [], 'Reside_time', [], 'GridInd_List', [], 'Arrival_Posx', [], 'Arrival_Posy', []);
    
    for i_part = 1: nparticles               
        Arrival_TF = (diff_GridInd(:, i_part) ~= 0);
        
        NArrival = sum(Arrival_TF);
        
        Arrival_time = zeros(NArrival+1, 1);
        Reside_time = zeros(NArrival+1, 1);
        GridInd_List = zeros(NArrival+1, 1);
        
        Arrival_Posx = zeros(NArrival+1, 1);
        Arrival_Posy = zeros(NArrival+1, 1);
        
        
        % absPos:  ||-t1-||-t2-||-t3-||-t4-|| (Resid_time = 4)
        % diffPos: XX-oo-u1-oo-u2-oo-u3-oo-XX
        
        % Initial position
        Arrival_time(1) = 1;
        GridInd_List(1) = GridInd(1, i_part);
        
        Arrival_Posx(1) = alphax(1, i_part);
        Arrival_Posy(1) = alphay(1, i_part);
        
        
        % Intermediate residence
        Reside_counter = 1;
        Arrival_counter = 1;
        for i_t = 2: Nts
          
            % Assistance:
            % x(i_t): i_t -> t_i
            % diffx(i_t) -> 1/2*((t_i) + (t_{i+1}));     (after t_i)
            % diffx(i_t-1) -> 1/2*((t_{i-1}) + (t_{i})); (before t_i)

            GridJump = diff_GridInd((i_t-1), i_part);
            if (GridJump ~= 0)
                Arrival_counter = Arrival_counter + 1;
                
                % Data-structure: List
                Arrival_time(Arrival_counter) = i_t;
                Reside_time(Arrival_counter-1) = Reside_counter;
                
                GridInd_List(Arrival_counter) = GridInd(i_t, i_part);
                
                Arrival_Posx(Arrival_counter) = alphax(i_t, i_part);
                Arrival_Posy(Arrival_counter) = alphay(i_t, i_part);
                
                Reside_counter = 1;
                
            else
                Reside_counter = Reside_counter + 1;
            end
        end
        assert(Arrival_counter == (NArrival+1));
        
        % Final residence
        Reside_time(end) = Reside_counter;
        
        % Storage
        Transit_Part(i_part).Arrival_time = Arrival_time;
        Transit_Part(i_part).Reside_time = Reside_time;
        Transit_Part(i_part).GridInd_List = GridInd_List;
        
        Transit_Part(i_part).Arrival_Posx = Arrival_Posx;
        Transit_Part(i_part).Arrival_Posy = Arrival_Posy;
    end
    
    
    % Summarise the particle-based data into cell-based data
    Transit_Stat = struct('Reside_time', [], ...
        'Reside_time_I', [], 'Reside_time_F', [], ...
        'GridInd_To_I', [],  'GridInd_From_F', [], ...
        'GridInd_From', [], 'GridInd_To', [], ...
        'Arrival_Posx', [], 'Arrival_Posy', [], ...
        'Exit_Posx', [], 'Exit_Posy', [], ...
        'Part_ID',[], 'Part_ID_F', []);
    
    Transit_Stat(length(Mesh)).GridInd_To = [];
    
    for i_part = 1: nparticles       
        Part_Reside_time = Transit_Part(i_part).Reside_time;
        Part_GridInd_List = Transit_Part(i_part).GridInd_List;
        Part_Arrival_Posx = Transit_Part(i_part).Arrival_Posx;
        Part_Arrival_Posy = Transit_Part(i_part).Arrival_Posy;
        
        % First Sample
        SampleInd = 1;
        
        cell_k = Part_GridInd_List(SampleInd);
        Transit_Stat(cell_k).Reside_time_I = [Transit_Stat(cell_k).Reside_time_I; Part_Reside_time(SampleInd)];
        Transit_Stat(cell_k).GridInd_To_I = [Transit_Stat(cell_k).GridInd_To_I; Part_GridInd_List(SampleInd+1)];
        
        
        % Exclude the final position
        %%% Can pre-calculate the size of array by looping over all particles
        for SampleInd = 2:(length(Part_GridInd_List)-1)
            cell_k = Part_GridInd_List(SampleInd);
            
            Transit_Stat(cell_k).Reside_time = [Transit_Stat(cell_k).Reside_time; Part_Reside_time(SampleInd)];
            Transit_Stat(cell_k).GridInd_From = [Transit_Stat(cell_k).GridInd_From; Part_GridInd_List(SampleInd-1)];
            Transit_Stat(cell_k).GridInd_To = [Transit_Stat(cell_k).GridInd_To; Part_GridInd_List(SampleInd+1)];
            
            Transit_Stat(cell_k).Arrival_Posx = [Transit_Stat(cell_k).Arrival_Posx; Part_Arrival_Posx(SampleInd)];
            Transit_Stat(cell_k).Arrival_Posy = [Transit_Stat(cell_k).Arrival_Posy; Part_Arrival_Posy(SampleInd)];
            
            Transit_Stat(cell_k).Exit_Posx = [Transit_Stat(cell_k).Exit_Posx; Part_Arrival_Posx(SampleInd+1)];
            Transit_Stat(cell_k).Exit_Posy = [Transit_Stat(cell_k).Exit_Posy; Part_Arrival_Posy(SampleInd+1)];
            
            Transit_Stat(cell_k).Part_ID = [Transit_Stat(cell_k).Part_ID; i_part];
        end
        
        % Final Cell
        SampleInd = length(Part_GridInd_List);
        
        cell_k = Part_GridInd_List(SampleInd);
        Transit_Stat(cell_k).Reside_time_F = [Transit_Stat(cell_k).Reside_time_F; Part_Reside_time(SampleInd)];
        Transit_Stat(cell_k).GridInd_From_F = [Transit_Stat(cell_k).GridInd_From_F; Part_GridInd_List(SampleInd-1)];
        
        Transit_Stat(cell_k).Part_ID_F = [Transit_Stat(cell_k).Part_ID_F; i_part];
    end
        
    save(filename, 'Analysed_Traj', 'Transit_Stat', 'Transit_Part', 'Mesh_Struct', '-v7.3')
    disp('Finished Saving Transition Statistics')
else
    load(filename)
end
%% Compute Auto-Correlation
Max_Reside_time = zeros(length(Mesh), 1);
Mode_Reside_time = zeros(length(Mesh), 1);
Median_Reside_time = zeros(length(Mesh), 1);
Mean_Reside_time = zeros(length(Mesh), 1);

% Displacement of particles
%assert(all(Analysed_Traj.x(:) == x(:)));

for cell_k = 1:length(Mesh)
    Transit_Stat_cell_k = Transit_Stat(cell_k);
    
    Max_Reside_time(cell_k) = max([Transit_Stat_cell_k.Reside_time(:); Transit_Stat_cell_k.Reside_time_F(:); Transit_Stat_cell_k.Reside_time_I(:)]);
    
    Mode_Reside_time(cell_k) = mode(Transit_Stat_cell_k.Reside_time(:));
    Median_Reside_time(cell_k) = median(Transit_Stat_cell_k.Reside_time(:));
    Mean_Reside_time(cell_k) = mean(Transit_Stat_cell_k.Reside_time(:));
end

ResideTime_Stat = struct('Max_Reside_time', Max_Reside_time, ...
                         'Mode_Reside_time', Mode_Reside_time, ...
                         'Median_Reside_time', Median_Reside_time, ...
                         'Mean_Reside_time', Mean_Reside_time);

Correleation_Stat = struct('Autocorr_SpatChunk', [], 'Autocorr_TempChunk', [], ...
    'Autocorr_SpatChunk_XX', [], 'Autocorr_TempChunk_XX', [], ...
    'Autocorr_SpatChunk_YY', [], 'Autocorr_TempChunk_YY', [], ...
    'Autocorr_SpatChunk_XY', [], 'Autocorr_TempChunk_XY', [], ...
    'AvailSample_SpatChunk', [],  'AvailSample_TempChunk', [] );
                     
TempChunk_HalfLength = 10;
for cell_k = 1:length(Mesh)
    Correleation_Stat(cell_k).Autocorr_SpatChunk = zeros(Max_Reside_time(cell_k)-1, 1);
    Correleation_Stat(cell_k).Autocorr_SpatChunk_XX = zeros(Max_Reside_time(cell_k)-1, 1);
    Correleation_Stat(cell_k).Autocorr_SpatChunk_YY = zeros(Max_Reside_time(cell_k)-1, 1);
    Correleation_Stat(cell_k).Autocorr_SpatChunk_XY = zeros(Max_Reside_time(cell_k)-1, 1);

    Correleation_Stat(cell_k).AvailSample_SpatChunk = zeros(Max_Reside_time(cell_k)-1, 1);
    
    
    Correleation_Stat(cell_k).Autocorr_TempChunk = zeros(TempChunk_HalfLength*2, 1);
    Correleation_Stat(cell_k).Autocorr_TempChunk_XX = zeros(TempChunk_HalfLength*2, 1);
    Correleation_Stat(cell_k).Autocorr_TempChunk_YY = zeros(TempChunk_HalfLength*2, 1);
    Correleation_Stat(cell_k).Autocorr_TempChunk_XY = zeros(TempChunk_HalfLength*2, 1);

    Correleation_Stat(cell_k).AvailSample_TempChunk = zeros(TempChunk_HalfLength*2, 1);
end

% Particle-base assemble of the correlation
Nts = length(Analysed_Traj.ts_list);
nparticles = size(Analysed_Traj.x, 2);

Part_GridInd = Analysed_Traj.GridInd;
x = Analysed_Traj.x;
y = Analysed_Traj.y;
diffx = diff(Analysed_Traj.x, 1);
diffy = diff(Analysed_Traj.y, 1);

for i_part = 1:nparticles
    % Explanation:
    % SpatChunk: only consider the trunk of sample trajectories in the bin
    % unbinned: consider fixed length sample trajectories
    
    Arrival_time = Transit_Part(i_part).Arrival_time;
    Reside_time = Transit_Part(i_part).Reside_time;

    % SpatChunk: consider trajectories fully in the bin
    for SpatChunk_ind = 1:length(Arrival_time)
        Chunk_starting_ind = Arrival_time(SpatChunk_ind);
        Chunk_ending_ind = Chunk_starting_ind + Reside_time(SpatChunk_ind)-1;  % In the time-series of position, not diff(position)!

        % During this interval the GridInd remains unchanged
        Chunk_absPos_ind = Chunk_starting_ind:Chunk_ending_ind;
        Chunk_GridInd = Part_GridInd(Chunk_starting_ind, i_part);
        
        assert(range(Part_GridInd(Chunk_absPos_ind, i_part))==0)
        if Chunk_ending_ind < Nts
            assert(Chunk_GridInd ~= Part_GridInd(Chunk_ending_ind+1, i_part))
        elseif Chunk_starting_ind > 1
            assert(Chunk_GridInd ~= Part_GridInd(Chunk_starting_ind-1, i_part))
        end

        Chunk_x = x(Chunk_absPos_ind, i_part); 
        Chunk_y = y(Chunk_absPos_ind, i_part);

        % Select the jumps between the chunk interval
        Chunk_diffPos_ind = Chunk_starting_ind:(Chunk_ending_ind-1);
        Chunk_diffx = diffx(Chunk_diffPos_ind, i_part);  % -1 again: to change from the view point of position to diff(position)
        Chunk_diffy = diffy(Chunk_diffPos_ind, i_part);  % -1 again: to change from the view point of position to diff(position)
        
        Chunk_diffr = sqrt(Chunk_diffx.^2 + Chunk_diffy.^2);
        
        % Sum the autocorrelation sample into stat
        Chunk_diffPos_length = length(Chunk_diffPos_ind);
        
        % absPos:  ||-t1-||-t2-||-t3-||-t4-|| (Resid_time = 4)
        % diffPos: XX-oo-u1-oo-u2-oo-u3-oo-XX
        for chunk_LocalInd = 1:Chunk_diffPos_length
            for tau = (0:(Chunk_diffPos_length - chunk_LocalInd))+1
                Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk(tau) + Chunk_diffr(chunk_LocalInd)*Chunk_diffr(chunk_LocalInd+(tau-1));
                
                Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk_XX(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk_XX(tau) + Chunk_diffx(chunk_LocalInd)*Chunk_diffx(chunk_LocalInd+tau-1);
                Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk_YY(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk_YY(tau) + Chunk_diffy(chunk_LocalInd)*Chunk_diffy(chunk_LocalInd+tau-1);
                Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk_XY(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_SpatChunk_XY(tau) + Chunk_diffx(chunk_LocalInd)*Chunk_diffy(chunk_LocalInd+tau-1);

                Correleation_Stat(Chunk_GridInd).AvailSample_SpatChunk(tau) = Correleation_Stat(Chunk_GridInd).AvailSample_SpatChunk(tau) + 1;
            end
        end
    end
    
    % TempChunk: consider trajectories fully in the bin
    for TempChunk_ind = 1:length(Arrival_time)
        ChunkBinReside_time = Reside_time(TempChunk_ind);
        Chunk_centre_ind = Arrival_time(TempChunk_ind) + floor(ChunkBinReside_time*0.5);
        
        Chunk_starting_ind = Chunk_centre_ind - TempChunk_HalfLength;
        Chunk_ending_ind = Chunk_centre_ind + TempChunk_HalfLength;  % In the time-series of position, not diff(position)!

        Chunk_starting_ind = max(Chunk_starting_ind, 1);  % Ensure not out of boundary
        Chunk_ending_ind = min(Chunk_ending_ind, Nts);
        
        % During this interval do not wish GridInd to vary too much...
        Chunk_absPos_ind = Chunk_starting_ind:Chunk_ending_ind;
        Chunk_GridInd = Part_GridInd(Chunk_starting_ind, i_part);

        Chunk_x = x(Chunk_absPos_ind, i_part); 
        Chunk_y = y(Chunk_absPos_ind, i_part);

        % Select the jumps between the chunk interval
        Chunk_diffPos_ind = Chunk_starting_ind:(Chunk_ending_ind-1);
        Chunk_diffx = diffx(Chunk_diffPos_ind, i_part);  % -1 again: to change from the view point of position to diff(position)
        Chunk_diffy = diffy(Chunk_diffPos_ind, i_part);  % -1 again: to change from the view point of position to diff(position)
        
        Chunk_diffr = sqrt(Chunk_diffx.^2 + Chunk_diffy.^2);
        
        % Sum the autocorrelation sample into stat
        Chunk_diffPos_length = length(Chunk_diffPos_ind);
        %assert(Chunk_diffPos_length ==  2*TempChunk_HalfLength);
        
        % absPos:  ||-t1-||-t2-||-t3-||-t4-|| (Resid_time = 4)
        % diffPos: XX-oo-u1-oo-u2-oo-u3-oo-XX
        for chunk_LocalInd = 1:Chunk_diffPos_length
            for tau = (0:(Chunk_diffPos_length - chunk_LocalInd))+1
                Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk(tau) + Chunk_diffr(chunk_LocalInd)*Chunk_diffr(chunk_LocalInd+(tau-1));
                
                Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk_XX(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk_XX(tau) + Chunk_diffx(chunk_LocalInd)*Chunk_diffx(chunk_LocalInd+tau-1);
                Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk_YY(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk_YY(tau) + Chunk_diffy(chunk_LocalInd)*Chunk_diffy(chunk_LocalInd+tau-1);
                Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk_XY(tau) = Correleation_Stat(Chunk_GridInd).Autocorr_TempChunk_XY(tau) + Chunk_diffx(chunk_LocalInd)*Chunk_diffy(chunk_LocalInd+tau-1);

                Correleation_Stat(Chunk_GridInd).AvailSample_TempChunk(tau) = Correleation_Stat(Chunk_GridInd).AvailSample_TempChunk(tau) + 1;
            end
        end
    end
    
    
    if mod(i_part, 100) == 0
        disp(['Finished particle ', num2str(i_part)]);
    end
end

DataTypeString = ['_CorreStat'];
filename = ['/home/s1046972/data/', TrajType_Save, DataTypeString, '_npart', num2str(nparticles), '_SamplingIntv', num2str(SamplingInterval_vis), '_Idx', num2str(Nx_cell) '.mat'];
    
save(filename, 'Correleation_Stat', 'ResideTime_Stat', '-v7.3')
disp('Finished Saving Autocorrelation Statistics')