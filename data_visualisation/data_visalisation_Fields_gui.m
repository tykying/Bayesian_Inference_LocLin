function data_visalisation_Fields_gui(arg)
global sliderVar popupVar
global sliderVar_cur

global ContPan
global figObj

addpath('../PlotUtilities/')
addpath('../Utilities/')
addpath('../StructConversion/')


%% List of Variables to Output 
% Explanation: 
% Var_Val: Values of Possible Parameters
% Var_name: Names of Possible Parameters
% Dynamic Variables will be created to based on 'Var_name'

% Parameters in sliders (Integer-valued, from 1 to max)
Layer_List = [2, 2];
Nx_cell_ARG_List = [8, 16, 32];
SamplingInterval_vis_List = [1, 2, 4, 8, 16, 32, 64, 128];
DS_rate_List = [4, 2, 1];
Nsteps_pc_List = [256, 512, 1024];
nTrial_List = [1:3];

% Parameters in popup menu
InfScheme_List = {'Local'};
veloc_Profile_List = {'QGM2_nReal_DSTemp', 'QGM2_nReal_DSSpat'};


%% Ad-hoc configuration to check sensitivity to initial conditions
% Parameters in sliders (Integer-valued, from 1 to max)
Layer_List = [2, 2];
Nx_cell_ARG_List = [16, 16];
SamplingInterval_vis_List = [1, 16:16:128, 160];
%SamplingInterval_vis_List = SIv_in;
DS_rate_List = [1,1];
Nsteps_pc_List = [8192*2, 8192*2];
nTrial_List = [1:3];

% Parameters in popup menu
InfScheme_List = {'Local'};
veloc_Profile_List = {'QGM2_DStemp'};


%% Fields to be visualised [Need to specially code-up what to visualise for each case]
Field_vis_List = {'BasicStat', 'Ellipse', 'Velocity', 'Diffusivity', ...
                    'Velocity_Init', 'Diffusivity_Init', ...
                    'Velocity_Final', 'Diffusivity_Final', ...
                    'Diffusivity_Fxy', 'TransDensity', ...
                    'theta', 'theta_Scale', ...
                    'PWC_Fields'};
plot_streamline = 0;
MH_step_List = [1:(32+1)];

% Lists = {Layer_List, Nx_cell_ARG_List, SamplingInterval_vis_List, DS_rate_List, Nsteps_pc_List, nTrial_List};
% for k = 1:length(Lists)
%     if length(Lists{k}) == 1
%         Lists{k} = [Lists{k}, Lists{k}];
%     end
% end

%% Set up sliders and popup menu
% Set up sliders
sliderVar_Range = {num2cell(Layer_List), num2cell(Nx_cell_ARG_List), ...
                num2cell(SamplingInterval_vis_List), ...
                num2cell(DS_rate_List), num2cell(Nsteps_pc_List), ...
                num2cell(nTrial_List), num2cell(MH_step_List)};
sliderVar_name = {'Layer', 'Nx_cell_ARG', 'SamplingInterval_vis', 'DS_rate', 'Nsteps_pc', 'nTrial', 'MH_step'};
assert(length(sliderVar_Range) == length(sliderVar_name));

% Set up popup menu
popupVar_Range = {InfScheme_List, veloc_Profile_List, Field_vis_List};
popupVar_name = {'InfScheme', 'veloc_Profile', 'Field_vis'};
assert(length(popupVar_Range) == length(popupVar_name));


%% Parameters for Graphics [No need to tune]
% Relative Positions
ContPan_dy = 0.15;

slider_xI = 0.25;
slider_dx = 0.20;
slider_dy = 1/length(sliderVar_name);
Vallabel_dx = 0.05;

slider_text_xI = slider_xI - 0.20;
textlabel_dx = 0.10;

curlabel_xI = slider_xI - textlabel_dx;

popup_xI = 0.75;
popup_dx = 0.20;
popup_dy = 1/length(popupVar_name);

popup_text_xI = popup_xI - 0.20;


% Specify the domain for plotting 
fig_xI = 0.05;
fig_yI = 0.15;
fig_dx = 1-fig_xI*2;
fig_dy = (1-ContPan_dy)-fig_yI*2;

% Specify the domain for explanation
exp_xI = 0.05;
exp_yI = 0.0;
exp_dy = fig_yI;


% If no argument input, set arg <- 1
if nargin == 0, arg = -1; end

switch arg
    
    case -1  %Initialisation
        warning off
        
        % Set seed to randn
        %randn('state',sum(100*clock))
        
        figObj = figure('Name','Posterior Mean','NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);
        
        % gcf = current figure handle
        %fig = gcf;
        clf
        % UserData: Memory space store in the object fig; enable communications
        set(gcf,'userdata',1);
        
        % Control Panel
        ContPan = uipanel(figObj, 'Title', 'Control Panel', ...
            'FontSize',12,...
            'Position',[0.0 1-ContPan_dy 1.0 ContPan_dy]);
        
        % Define sliders (ONLY ALLOW INTEGER VALUE FOR VISUALISATION)
        % Value: initial value
        % CallBack: whenever the slider is altered, call function 'data_visalisation_Fields_gui(1)'
        for k = 1:length(sliderVar_name)
            bar_yI = 1-k*slider_dy;
            sliderVar{k} = uicontrol('Parent',ContPan, 'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [slider_xI bar_yI slider_dx slider_dy], ...
                'Min', 1, 'Max', length(sliderVar_Range{k}), ...
                'Value', 1, ...
                'CallBack', 'data_visalisation_Fields_gui(1)', ...
                'SliderStep', [1/(length(sliderVar_Range{k})-1) 1]);
            
            % Define the text controls for the minimum and maximum values
%             % Draw a text at the designated position
%             sliderVar_min{k} = uicontrol('Parent',ContPan, 'Style','text', ...
%                 'Units', 'normalized', ...
%                 'Pos', [slider_xI-Vallabel_dx bar_yI Vallabel_dx slider_dy], ...
%                 'String',num2str(get(sliderVar{k}, 'Min')));
%             
%             sliderVar_max{k} = uicontrol('Parent',ContPan, 'Style', 'text', ...
%                 'Units', 'normalized', ...
%                 'Pos', [slider_xI+slider_dx bar_yI Vallabel_dx slider_dy], ...
%                 'String',num2str(get(sliderVar{k}, 'Max')));
            
            
            % Define the slider labels
            sliderVar_label{k} = uicontrol('Parent',ContPan, 'Style', 'text', ...
                'Units','normalized', ...
                'Pos', [slider_text_xI bar_yI textlabel_dx slider_dy], ...
                'String', sliderVar_name{k});
            
            % Define the text controls for the current values
            sliderVar_cur{k} = uicontrol('Parent',ContPan, 'Style', 'text', ...
                'Units', 'normalized', ...
                'Pos', [curlabel_xI bar_yI Vallabel_dx slider_dy],...
                'String',num2str(sliderVar_Range{k}{get(sliderVar{k},'Val')}));
            
        end
        
        for k = 1:length(popupVar_name)
            popup_yI = 1-k*popup_dy;
            popupVar{k} = uicontrol('Parent',ContPan, 'Style', 'popupmenu', ...
                'Units', 'normalized', ...
                'Position', [popup_xI popup_yI popup_dx popup_dy], ...
                'String', popupVar_Range{k}, ...
                'CallBack', 'data_visalisation_Fields_gui(1)');
            
            % Define the slider labels
            popupVar_label{k} = uicontrol('Parent',ContPan, 'Style', 'text', ...
                'Units','normalized', ...
                'Pos', [popup_text_xI popup_yI textlabel_dx popup_dy], ...
                'String', popupVar_name{k});
        end
        
        % Draw now
        data_visalisation_Fields_gui(0);
        
%         % Define push buttons
%         
%         pb_info = uicontrol(fig,'Style','push','Units','normalized',...
%             'Position',[0.85 0.43 0.125 0.07],'String','Info',...
%             'Callback','helpwin data_visalisation_Fields_gui');
%         
%         pn_clw = uicontrol(gcf,'Style','push', 'Units','normalized',...
%             'Pos',[0.85 0.34 0.125 0.07],'String','Close','Callback','close(gcf),warning on');   
         
    case 0 % animate
        % Dynamic Variables will be created to based on 'Var_name'
        Var_List = [sliderVar, popupVar];
        VarName_Range = [sliderVar_Range, popupVar_Range];  % _Val: Includes all possible values
        VarName_List = [sliderVar_name, popupVar_name];
        run('../Scripts/Script_VarName_List_to_Var');
        
        
        % If want to use for other data: edit from here             
        
        if (isempty(strfind(Field_vis, 'BasicStat'))==0)  % =0: Matched
            DataTypeString = 'BinnedTrajData';
            run('../Scripts/Script_Filenames');
            BinnedTrajData = load(filename_BTD);
            
            DataTypeString = 'BayesSampleStat';
            run('../Scripts/Script_Filenames');
            BayesSampleStat = load(filename_BSS);
        
            RectMesh_Param = BinnedTrajData.Mesh_Struct.RectMesh_Param;
            Mesh = BinnedTrajData.Mesh_Struct.Mesh;
            
            JumpsCFL = BinnedTrajData.TrajJumps_DA.JumpsCFL;
            NJumps_cell = BayesSampleStat.MCMC_Param.NJumps_cell;
            
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
            NJumps_cell_k = convert_cellk_to_cellij( NJumps_cell, RectMesh_Param );
            JumpsCFL_k = convert_cellk_to_cellij( JumpsCFL, RectMesh_Param );
            
            contourf_vis_data = { JumpsCFL_k(:,:,1), JumpsCFL_k(:,:,2); NJumps_cell_k, max(JumpsCFL_k(:,:,1), JumpsCFL_k(:,:,2))};
            
            title_data = {'JumpsCFL: x', 'JumpsCFL: y'; 'NJumps', 'JumpsCFL: max(x,y)'};
            colormap_data = {'jet', 'jet'; 'jet', 'jet'};
            caxis_data = {[0, 1], [0, 1]; [], [0, 1.25]};

            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
            
            ax = pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);
        end
        
        
        if (isempty(strfind(Field_vis, 'theta'))==0)  % =0: Matched
            DataTypeString = 'BayesSampleStat';
            run('../Scripts/Script_Filenames');
            BayesSampleStat = load(filename_BSS);
        
            RectMesh_Param = BayesSampleStat.Mesh_Struct.RectMesh_Param;
            Mesh = BayesSampleStat.Mesh_Struct.Mesh;
            
            theta_Stat = BayesSampleStat.theta_Stat;
            
            theta_store = theta_Stat.theta_store;
                       
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
            
            step = MH_step;
            theta_step = squeeze(theta_store(:,:,step));
            theta_step_k = convert_cellk_to_cellij( theta_step, RectMesh_Param );

            % u, v, Sigma1, Sigma 2 at cell centres
            contourf_vis_data = {theta_step_k(:,:,1), theta_step_k(:,:,2), ...
                theta_step_k(:,:,7), theta_step_k(:,:,8);  ...
                theta_step_k(:,:,3), theta_step_k(:,:,4), ...
                theta_step_k(:,:,5), theta_step_k(:,:,9); };
            
            % Set visualisation
            if (isempty(strfind(veloc_Profile, 'QGM2')) == 0)
                UV_lim = 0.1;
                sigma_lim = 125;
                
                dudx_range = [-1.5E-6 1.5E-6];
                dudy_range = [-1.5E-6 1.5E-6];
                dvdx_range = [-1.5E-6 1.5E-6];
            else
                UV_lim = 1.5;
                sigma_lim = 0.75;
            end
            
            U_range = [-UV_lim/5*8 UV_lim/5*8];
            U_CLV = linspace(-UV_lim, UV_lim, 21);           
            
            sigma_range = [0 sigma_lim];
            phi_range = [-pi/2 pi/2];
            sigma_CLV = linspace(0, sigma_lim, 21);
            phi_CLV = [-pi/2: pi/18: pi/2];
            
            title_data = {'u', 'v', '\sigma_1', '\sigma_2'; 'du/dx', 'du/dy', 'dv/dx', '\phi'};
            colormap_data = {'jet', 'jet', 'jet', 'jet'; 'jet', 'jet', 'jet', 'hsv' };
            caxis_data = {U_range, U_range, sigma_range, sigma_range;
                          dudx_range, dudy_range, dvdx_range, phi_range};
                      
            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
            
            ax = pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);
        end
        
        if (isempty(strfind(Field_vis, 'theta_Scale'))==0)  % =0: Matched
            DataTypeString = 'BayesSampleStat';
            run('../Scripts/Script_Filenames');
            BayesSampleStat = load(filename_BSS);
        
            RectMesh_Param = BayesSampleStat.Mesh_Struct.RectMesh_Param;
            Mesh = BayesSampleStat.Mesh_Struct.Mesh;
            
            theta_Scale = BayesSampleStat.theta_Stat.theta_Scale;
            
                       
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
            
            theta_Scale_k = convert_cellk_to_cellij( theta_Scale, RectMesh_Param );

            % u, v, Sigma1, Sigma 2 at cell centres
            contourf_vis_data = {theta_Scale_k(:,:,1), theta_Scale_k(:,:,2), ...
                theta_Scale_k(:,:,7), theta_Scale_k(:,:,8);  ...
                theta_Scale_k(:,:,3), theta_Scale_k(:,:,4), ...
                theta_Scale_k(:,:,5), theta_Scale_k(:,:,9); };
            
            % Set visualisation
            if (isempty(strfind(veloc_Profile, 'QGM2')) == 0)
                UV_lim = 0.1;
                sigma_lim = 125;
                
                dudx_range = [0 max(max(theta_Scale_k(:,:,3)))];
                dudy_range = [0 max(max(theta_Scale_k(:,:,4)))];
                dvdx_range = [0 max(max(theta_Scale_k(:,:,5)))];
            else
                UV_lim = 1.5;
                sigma_lim = 0.75;
            end
            
            U_range = [-UV_lim/5*8 UV_lim/5*8];
            U_CLV = linspace(-UV_lim, UV_lim, 21);           
            
            sigma_range = [0 sigma_lim];
            phi_range = [-pi/2 pi/2];
            sigma_CLV = linspace(0, sigma_lim, 21);
            phi_CLV = [-pi/2: pi/18: pi/2];
            
            title_data = {'u', 'v', '\sigma_1', '\sigma_2'; 'du/dx', 'du/dy', 'dv/dx', '\phi'};
            colormap_data = {'jet', 'jet', 'jet', 'jet'; 'jet', 'jet', 'jet', 'hsv' };
            caxis_data = {U_range, U_range, sigma_range, sigma_range;
                          dudx_range, dudy_range, dvdx_range, phi_range};
                      
            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
            
            ax = pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);
        end
        
        if  (isempty(strfind(Field_vis, 'Velocity'))==0)  % =0: Matched
            DataTypeString = 'BayesSampleMeanFine';
            run('../Scripts/Script_Filenames');
            BayesSampleMeanFine = load(filename_BSMF);
            
            grid_vis_f = BayesSampleMeanFine.grid_vis_f;
            
            if (isempty(strfind(Field_vis, '_Init'))==0)
                U_vis = BayesSampleMeanFine.theta_Init_f.u;
                V_vis = BayesSampleMeanFine.theta_Init_f.v;
            elseif (isempty(strfind(Field_vis, '_Final'))==0)
                U_vis = BayesSampleMeanFine.theta_Final_f.u;
                V_vis = BayesSampleMeanFine.theta_Final_f.v;
            else
                U_vis = BayesSampleMeanFine.MeanFields_f.u;
                V_vis = BayesSampleMeanFine.MeanFields_f.v;
            end
            
            % Set visualisation
            if (isempty(strfind(veloc_Profile, 'QGM2')) == 0)
                UV_lim = 0.3;
            else
                UV_lim = 1.5;
            end
            
            U_range = [-UV_lim/5*8 UV_lim/5*8];
            U_CLV = linspace(-UV_lim, UV_lim, 21);            
            
            contourf_vis_data = {U_vis, V_vis};
            ContourLineValue_data = {U_CLV, U_CLV};
            caxis_data = {U_range, U_range};
            colormap_data = {'jet', 'jet'};
            
            title_data = {'u', 'v'};
            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'ContourLineValue_data', ContourLineValue_data, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
                        
            ax = contourf_subplot(grid_vis_f, contourf_vis_data, subplot_param);
            
            % Add streamline on top of the velocity field
            if plot_streamline == 1
                nPlot_y = size(ax, 1);
                nPlot_x = size(ax, 2);
                
                for j = 1:nPlot_x
                    for i = 1:nPlot_y
                        plot_MeanFlowStreamline(ax{i, j}, grid_vis_f, U_vis, V_vis);
                    end
                end
            end
            
        end
        
        if (isempty(strfind(Field_vis, 'Diffusivity'))==0)  % =0: Matched
            DataTypeString = 'BayesSampleMeanFine';
            run('../Scripts/Script_Filenames');
            BayesSampleMeanFine = load(filename_BSMF);
            
            grid_vis_f = BayesSampleMeanFine.grid_vis_f;
            
            if (isempty(strfind(Field_vis, '_Init'))==0)
                K_sigma1_vis = BayesSampleMeanFine.theta_Init_f.K_sigma1;
                K_sigma2_vis = BayesSampleMeanFine.theta_Init_f.K_sigma2;
                K_phi_vis = BayesSampleMeanFine.theta_Init_f.K_phi;
            elseif (isempty(strfind(Field_vis, '_Final'))==0)
                K_sigma1_vis = BayesSampleMeanFine.theta_Final_f.K_sigma1;
                K_sigma2_vis = BayesSampleMeanFine.theta_Final_f.K_sigma2;
                K_phi_vis = BayesSampleMeanFine.theta_Final_f.K_phi;
            elseif (isempty(strfind(Field_vis, '_Fxy'))==0)
                K_sigma1 = BayesSampleMeanFine.theta_Final_f.K_sigma1;
                K_sigma2 = BayesSampleMeanFine.theta_Final_f.K_sigma2;
                K_phi_vis = BayesSampleMeanFine.theta_Final_f.K_phi;
                
                [Kxx, Kyy, Kxy]  = Kpolar_to_Kcart_vectorised(K_sigma1, K_sigma2, K_phi_vis);
                K_sigma1_vis = sqrt(Kxx);
                K_sigma2_vis = sqrt(Kyy);
            else
                K_sigma1_vis = BayesSampleMeanFine.MeanFields_f.K_sigma1;
                K_sigma2_vis = BayesSampleMeanFine.MeanFields_f.K_sigma2;
                K_phi_vis = BayesSampleMeanFine.MeanFields_f.K_phi;
            end
            
            % Set visualisation
            if (isempty(strfind(veloc_Profile, 'QGM2')) == 0)
                sigma_lim = 125;                
            else
                sigma_lim = 1.0;
            end           
            sigma_range = [0 sigma_lim];
            phi_range = [-pi/2 pi/2];
            sigma_CLV = linspace(0, sigma_lim, 21);
            phi_CLV = [-pi/2: pi/18: pi/2];
                            
            contourf_vis_data = {K_sigma1_vis, K_sigma2_vis, K_phi_vis};
            ContourLineValue_data = {sigma_CLV, sigma_CLV, phi_CLV};
            caxis_data = {sigma_range, sigma_range, phi_range};
            colormap_data = {'jet', 'jet', 'hsv'};
            
            if (isempty(strfind(Field_vis, '_Fxy'))==0)
                title_data = {'\sqrt{K_{xx}}', '\sqrt{K_{yy}}', '\phi'};
            else
                title_data = {'\sigma_1', '\sigma_2', '\phi'};
            end


            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'ContourLineValue_data', ContourLineValue_data, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
                        
            contourf_subplot(grid_vis_f, contourf_vis_data, subplot_param);
        end
        
        if  (isempty(strfind(Field_vis, 'TransDensity_BACKUP'))==0)  % =0: Matched
            DataTypeString = 'BinnedTrajTransMat';
            run('../Scripts/Script_Filenames');
            BinnedTrajTransMat = load(filename_BTTM);
            
            grid_vis_f = BinnedTrajTransMat.grid_vis_f;
            
            pi_Stat_f = BinnedTrajTransMat.pi_Stat_f;
            RatioRemain_f = BinnedTrajTransMat.RatioRemain_f;
            RatioRemainNeigh_f = BinnedTrajTransMat.RatioRemainNeigh_f;
            
            Ncells = length(BinnedTrajTransMat.Mesh_Struct.Mesh);
            NJumps_pc = sum(sum(BinnedTrajTransMat.TransMat, 2))/Ncells;
            
            % Deviation from the expected stationary distribution
            Stationary_Distri_Anomaly = 1 - Ncells*pi_Stat_f;
            
            StatDist_Anomaly_range = [-0.15, 0.15];
            StatDist_Anomaly_CLV = [-0.15:0.015:0.15];
            Prob_range = [0 1];
            Prob_CLV = [0:0.1:1];
            
            contourf_vis_data = {Stationary_Distri_Anomaly, RatioRemain_f, RatioRemain_f+RatioRemainNeigh_f};
            ContourLineValue_data = {StatDist_Anomaly_CLV, Prob_CLV, Prob_CLV};
            caxis_data = {StatDist_Anomaly_range, Prob_range, Prob_range};
            colormap_data = {'jet', 'jet', 'jet'};
            
            title_data = {'Stationary Distribution', {['Jumps/Cell: ', num2str(NJumps_pc)]; 'Prob(Remain)'}, 'Prob(Remain+Neighbour)'};
            
            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'ContourLineValue_data', ContourLineValue_data, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
                        
            contourf_subplot(grid_vis_f, contourf_vis_data, subplot_param);
        end
        
        if  (isempty(strfind(Field_vis, 'TransDensity'))==0)  % =0: Matched
            DataTypeString = 'BinnedTrajTransMat';
            run('../Scripts/Script_Filenames');
            BinnedTrajTransMat = load(filename_BTTM);
            
            Mesh = BinnedTrajTransMat.Mesh_Struct.Mesh;
            RectMesh_Param = BinnedTrajTransMat.Mesh_Struct.RectMesh_Param;
            
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);

            pi_Stat = convert_cellk_to_cellij( BinnedTrajTransMat.pi_Stat', RectMesh_Param );
            RatioRemain = convert_cellk_to_cellij( BinnedTrajTransMat.RatioRemain, RectMesh_Param );
            RatioRemainNeigh = convert_cellk_to_cellij( BinnedTrajTransMat.RatioRemainNeigh, RectMesh_Param );
            
            Ncells = length(BinnedTrajTransMat.Mesh_Struct.Mesh);
            NJumps_pc = sum(sum(BinnedTrajTransMat.TransMat, 2))/Ncells;
            
            % Deviation from the expected stationary distribution
            Stationary_Distri_Anomaly = 1 - Ncells*pi_Stat;
            
            StatDist_Anomaly_range = [-0.15, 0.15];
            StatDist_Anomaly_CLV = [-0.15:0.015:0.15];
            Prob_range = [0 1];
            Prob_CLV = [0:0.1:1];
            
            contourf_vis_data = {Stationary_Distri_Anomaly, RatioRemain, RatioRemain+RatioRemainNeigh};
            ContourLineValue_data = {StatDist_Anomaly_CLV, Prob_CLV, Prob_CLV};
            caxis_data = {StatDist_Anomaly_range, Prob_range, Prob_range};
            colormap_data = {'jet', 'jet', 'jet'};
            
            title_data = {'Stationary Distribution', {['Jumps/Cell: ', num2str(NJumps_pc)]; 'Prob(Remain)'}, 'Prob(Remain+Neighbour)'};
            
            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'ContourLineValue_data', ContourLineValue_data, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
                        
            pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);
            
        end
        
        
        if (isempty(strfind(Field_vis, 'PWC_Fields'))==0)  % =0: Matched
            DataTypeString = 'BayesSampleStat';
            run('../Scripts/Script_Filenames');
            BayesSampleStat = load(filename_BSS);
            theta_Stat = BayesSampleStat.theta_Stat;
            theta_final = theta_Stat.theta_store(:, :, end);
            
            kappa_final = 0.5*(theta_final(:, 7).^2 + theta_final(:, 8).^2);
            
            DataTypeString = 'BinnedTrajData';
            run('../Scripts/Script_Filenames');
            BinnedTrajData = load(filename_BTD);
                        
            RectMesh_Param = BinnedTrajData.Mesh_Struct.RectMesh_Param;
            Mesh = BinnedTrajData.Mesh_Struct.Mesh;
            TrajJumps_MomentGlobal = BinnedTrajData.TrajJumps_DA.TrajJumps_MomentGlobal;
                        
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
            
            % Crude Estimates of Velocities: by simple mean of jumps; 
            % m1, m2 = 1st and 2nd moment
            U_m1jumps = TrajJumps_MomentGlobal.diffx(:, 1)./TrajJumps_MomentGlobal.h(:,1);
            V_m1jumps = TrajJumps_MomentGlobal.diffy(:, 1)./TrajJumps_MomentGlobal.h(:,1);
            
            U_m2jumps = TrajJumps_MomentGlobal.diffx(:, 2);
            V_m2jumps = TrajJumps_MomentGlobal.diffy(:, 2);
            
            kappa = 0.5*(U_m2jumps + V_m2jumps)./(2*TrajJumps_MomentGlobal.h(:,1));
            
            U_m1jumps_k = convert_cellk_to_cellij( U_m1jumps, RectMesh_Param );
            V_m1jumps_k = convert_cellk_to_cellij( V_m1jumps, RectMesh_Param );
            kappa_k = convert_cellk_to_cellij( kappa, RectMesh_Param );
            kappa_final_k = convert_cellk_to_cellij( kappa_final, RectMesh_Param );
            
            
            contourf_vis_data = {U_m1jumps_k, V_m1jumps_k, kappa_k, kappa_k-kappa_final_k};
            
            title_data = {'U', 'V', 'kappa', 'kappa pwl'};
            colormap_data = {'jet', 'jet', 'jet', 'jet'};
            caxis_data = {[], [], [0, 10000], [-250, 250]};

            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
            
            ax = pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);
        end

                
        if (isempty(strfind(Field_vis, 'Ellipse'))==0)  % =0: Matched
            DataTypeString = 'BayesSampleStat';
            run('../Scripts/Script_Filenames');
            BayesSampleStat = load(filename_BSS);
        
            RectMesh_Param = BayesSampleStat.Mesh_Struct.RectMesh_Param;
            Mesh = BayesSampleStat.Mesh_Struct.Mesh;
            
            theta_Stat = BayesSampleStat.theta_Stat;
            MCMC_NSample = BayesSampleStat.MCMC_Stat.MCMC_NSample;
            theta_mean = theta_Stat.theta_sum/MCMC_NSample;
            thetasq_mean = theta_Stat.thetasq_sum/MCMC_NSample;
                       
            grid_vis = setup_grid_vis(Mesh, RectMesh_Param);            
            
            u = theta_mean(:, 1);
            v = theta_mean(:, 2);
            sigma1 = theta_mean(:, 7);
            sigma2 = theta_mean(:, 8);
            kappa = 0.5*(sigma1.^2 + sigma2.^2);
            phi = theta_mean(:, 9);
            
            eccentricity2 = sqrt((sigma1./sigma2).^2 - 1);
            
            kappa_ij = convert_cellk_to_cellij( kappa, RectMesh_Param );
            u_ij = convert_cellk_to_cellij( u, RectMesh_Param );
            v_ij = convert_cellk_to_cellij( v, RectMesh_Param );
            
            % Set visualisation
            if (isempty(strfind(veloc_Profile, 'QGM2')) == 0)
                sigma_lim = 100^2;                
            else
                sigma_lim = 1.0;
            end           
            sigma_range = [0 sigma_lim];
            phi_range = [-pi/2 pi/2];
            sigma_CLV = linspace(0, sigma_lim, 21);
            phi_CLV = [-pi/2: pi/18: pi/2];

            
            eccentricity2_ij = convert_cellk_to_cellij( eccentricity2, RectMesh_Param );
            ecc2_range = [0, 2.5];
            ecc2_CLV = linspace(0, 2.5, 21);

            contourf_vis_data = {kappa_ij, eccentricity2_ij, zeros(size(eccentricity2_ij))};
            ContourLineValue_data = {sigma_CLV, ecc2_CLV, []};
            caxis_data = {sigma_range, ecc2_range, [0,1]};
            colormap_data = {'jet','copper', 'gray'};
            
            title_data = {'1/2tr(kappa)', 'eccentricity', 'streamline'};
            
            figPos = [fig_xI fig_yI fig_dx fig_dy];
            subplot_param = struct('Position', figPos, ...
                'ContourLineValue_data', ContourLineValue_data, ...
                'caxis_data', caxis_data, 'title_data', title_data, ...
                'colormap_data', colormap_data);
            
            ax_subplot = pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);
                        
            %hold(ax_subplot{1,1}, 'on');
            plot_MeanFlowStreamline(ax_subplot{1, 3}, grid_vis, u_ij, v_ij);
            %quiver(ax_subplot{1,1}, grid_vis.X, grid_vis.Y, u_ij, v_ij)
            
            % Plot Ellipse
            dx_max = max(max(diff(grid_vis.cell_centre(:,:,1), 1)));
            dy_max = max(max(diff(grid_vis.cell_centre(:,:,2)', 1)));
            
            for cell_k = 1:length(Mesh)
                hold(ax_subplot{1,3}, 'on');
                
                sigma1_cellk = sigma1(cell_k);
                sigma2_cellk = sigma2(cell_k);
                phi_cellk = phi(cell_k);
                s = sqrt(dx_max*dy_max/(2*pi*sigma1_cellk*sigma2_cellk));
                a = s*sigma1_cellk;
                b = s*sigma2_cellk;
                centre = Mesh(cell_k).cell_centre;            
                PlotSpec = struct('Color', 'yellow', 'LineStyle', '-');
                ellipse_Param = struct('a', a, 'b', b, 'phi', phi_cellk, ...
                                    'centre', centre, 'PlotSpec', PlotSpec);
                PlotEllipse(ax_subplot{1,3}, ellipse_Param);
            end

        end
        
        
    case 1   % Alter slider values
        for k = 1:length(sliderVar)
            set(sliderVar_cur{k}, 'String', num2str(sliderVar_Range{k}{get(sliderVar{k},'Val')}) );
        end
        
        % Re-Draw
        data_visalisation_Fields_gui(0);
end

end
