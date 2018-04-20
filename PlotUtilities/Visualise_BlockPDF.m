%% Visualisation of data in blocks
function [SampleDataStat] = Visualise_BlockPDF(SampleData, SampleDataLabel, Block_Param, RectMesh_Param, PlotString_Param)

% Number of cells in a block
BlockLength_x = Block_Param.BlockLength_x;
BlockLength_y = Block_Param.BlockLength_y;
fit_curve = Block_Param.fit_curve;

Nx_cell = RectMesh_Param.Nx_cell;
Ny_cell = RectMesh_Param.Ny_cell;

blocksize = BlockLength_x*BlockLength_y;

% Number of block along each direction
NBlock_x = Nx_cell/BlockLength_x;
NBlock_y = Ny_cell/BlockLength_y;

NBlock = NBlock_x*NBlock_y;

SampleDataFields = fieldnames(SampleData(1));
N_Fields = numel(SampleDataFields);

N_Fields_vis = numel(SampleDataLabel);

N_Samples = numel(SampleData);
Ncell_SpDis = size(SampleData(1).(SampleDataFields{1}), 1);

assert(N_Fields >= N_Fields_vis);
assert(Ncell_SpDis == Nx_cell*Ny_cell);

% Copy the Parameters as local variables
Trackscheme = PlotString_Param.Trackscheme;
DGscheme = PlotString_Param.DGscheme;
NJumps_vis = PlotString_Param.NJumps_vis;
SamplingInterval_vis = PlotString_Param.SamplingInterval_vis;
FigTitle_suffix_line1 = PlotString_Param.FigTitle_suffix_line1;
FigTitle_suffix_line2 = PlotString_Param.FigTitle_suffix_line2;
traj_outputfolder = PlotString_Param.traj_outputfolder;
block_outputfolder = PlotString_Param.block_outputfolder;
fieldname = PlotString_Param.fieldname;
figure_vis = PlotString_Param.figure_vis;

% Compute SampleDataStat (store same way as theta_store)
SampleDataStat = Compute_SampleDataStat(SampleData);


% Check Validity of SampleDataLabel
N_ValidLabel = 0;
for j = 1: N_Fields
    FldName = SampleDataFields{j};
    for i = 1: N_Fields_vis
        LblName = SampleDataLabel{i};
        
        if strcmp(FldName, LblName)
            N_ValidLabel = N_ValidLabel + 1;
        end
    end
end
assert(N_ValidLabel == numel(SampleDataLabel))

% Visualisation in blocks
for block_j = 1:NBlock_y
    for block_i = 1:NBlock_x
        block_k = (block_j-1)*NBlock_x + block_i;
        
        % Calculate the bottomleft cell index in block k
        cell_BlockBegin_i = (block_i-1)*BlockLength_x + 1;
        cell_BlockBegin_j = (block_j-1)*BlockLength_y + 1;
        
        f100 = figure('visible', figure_vis);
        set(f100, 'Position', [100, 100, 800, 800]);
        
        subplot_ind = 0;
        
        % Reverse the direction of looping over cells in y-direction within
        % block: to make subplot align with the physical coordinates
        for cellj_BlockLoc = BlockLength_y:-1:1
            for celli_BlockLoc = 1:BlockLength_x
                cell_i = cell_BlockBegin_i + (celli_BlockLoc-1);
                cell_j = cell_BlockBegin_j + (cellj_BlockLoc-1);
                
                % Finally obtain the global index for the target cell
                cell_k = RectMesh_Param.cellij_to_k(cell_i, cell_j);
                
                subplot_ind = subplot_ind + 1;
                
                subplot(BlockLength_x, BlockLength_y, subplot_ind)
                % %                 colorbar
                % %                 xlabel('x (in km)')
                % %                 ylabel('y (in km)')
                % %                 %pbaspect([1 1 1])
                %                 axis off
                %                 marker = '.';
                %                 scatter(CellsJumps(cell).y, CellsJumps(cell).v, marker)
                %                 hold on
                %                 marker = 'o';
                %                 colour = 'red';
                %                 scatter(mean(CellsJumps(cell).y), mean(CellsJumps(cell).v), marker, colour)
               
                % A cell of histograms
                CellHistgram = cell(N_Fields_vis, 2);
                
                for i = 1: N_Fields_vis
                    SampleDataField = SampleDataLabel{i};
                    % Obtain the correct data structure for the Data
                    [RawDataNCell, RawDataDim] = size(SampleData(1).(SampleDataField));  % First dimension = number of cells
                    
                    CellRawData = zeros(N_Samples, RawDataDim);
                    for Sample = 1:N_Samples
                        CellRawData(Sample, :) = SampleData(Sample).(SampleDataField)(cell_k, :);
                    end
                    
                    for compon = 1:RawDataDim
                        nbins = 20;
                        CellHistgram{i, compon} = histogram(CellRawData(:,compon), nbins);
                        CellHistgram{i, compon}.Normalization = 'pdf';
                        hold on
                        
                        
                        CellRawData_Mean = SampleDataStat.(SampleDataField).Mean(cell_k,compon);
                        CellRawData_Vari = SampleDataStat.(SampleDataField).Vari(cell_k,compon);
                        CellRawData_Skew = SampleDataStat.(SampleDataField).Skew(cell_k,compon);
                        CellRawData_Kurt = SampleDataStat.(SampleDataField).Kurt(cell_k,compon);
                        
                        if fit_curve == 1
                             L_gau = ceil(max(1.2*abs(CellRawData(:, compon))));
                             
                             x_gau = linspace(-L_gau, L_gau, 1000);
                             norm = normpdf(x_gau, CellRawData_Mean, sqrt(CellRawData_Vari));
                             plot(x_gau, norm, 'red')
                             
                             fitted_Param = [CellRawData_Mean, CellRawData_Vari];
                             
                             title_string = ['Normal Fit: ', sprintf('(%4.2f, %4.2f)', fitted_Param(1), fitted_Param(2))];
                             title(title_string);
                        end

                        if fit_curve == 2
                            L_iga = max(1.2*abs(CellRawData(:, compon)));
                            x_iga = linspace(0, L_iga, 1000);
                            alpha = 1/(CellRawData_Vari/(CellRawData_Mean^2))+2;
                            beta = CellRawData_Mean * (alpha-1);
                            [ iga ] = inversegampdf( x_iga, alpha, beta );
                            plot(x_iga, iga, 'red')
                            
                            fitted_Param = [alpha, beta];
                            
                            title_string = ['IG Fit: ', sprintf('(%4.1f, %4.1f)', fitted_Param(1), fitted_Param(2))];
                            title(title_string);
                        end
                        
                        
                    end
                end
                %legend(SampleDataLabel)
                %legend('Sigma 1')
            end
        end
                
        FigTitle = {['P.D.F.', FigTitle_suffix_line1]; ...
            [FigTitle_suffix_line2, ...
            ['; Block i = ', num2str(block_i),'; Block j = ', num2str(block_j)]]};
        suptitle(FigTitle)
        
        % Configuration for output to eps
        param_list = ['_', Trackscheme, '_h', num2str(SamplingInterval_vis), 'hr', '_NJ' , num2str(NJumps_vis), 'K', '_', DGscheme, '_Idx', num2str(Nx_cell), '_Blc', num2str(block_k)];
        
        traj_outputpath = [block_outputfolder, 'pdf_', fieldname, param_list ,'.eps'];
        saveas(f100, traj_outputpath,'eps');
        
    end
end

end