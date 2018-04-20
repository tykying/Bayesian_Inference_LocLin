%% For Picked Theta
addpath('./PlotUtilities/')
grid_vis = setup_grid_vis(Mesh, RectMesh_Param);

% For Picked Theta
figure(400)
theta_Init_Picked_ij = convert_cellk_to_cellij( theta_Init_Picked, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_Init_Picked_ij);

figure(500)
theta_SSD_Picked_ij = convert_cellk_to_cellij( theta_SSD_Picked, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, theta_SSD_Picked_ij);

figure(600)
AccRatio_Picked_ij = convert_cellk_to_cellij( AccRatio_Picked, RectMesh_Param );
ax_sp = plot_dataij_k(grid_vis, AccRatio_Picked_ij);

%%

%    theta_Picked = zeros( size(theta_Init_List, 1), size(theta_Init_List, 2) );
%     AccRatio_Pick = zeros(size(theta_Picked));
%     nTrial_Picked = zeros(size(theta_Picked));
%     theta_SSD_Picked = zeros(size(theta_Picked));
%     
%     
%     for cell_k = 1:size(theta_Init_List, 1)
% %     end
% %     for cell_k = cell_k_List_TwoDiverged'
% 
%             theta_store_Mean_List_cell_k = squeeze(theta_Init_List(cell_k,:,:));
%             theta_SSD_List_cell_k = squeeze(theta_SSD_List(cell_k,:,:));
%             AccRatio_List_cell_k = squeeze(AccRatio_List(cell_k,:,:));
% 
%             % 0.23: benchmark
%             % Maximum deviation from the benchmark among all components
%             AccRatioRel_max_over_AllComp = max(abs(AccRatio_List_cell_k - 0.23), [], 1);
%             nTrial_maxAccRat = find(AccRatioRel_max_over_AllComp == min(AccRatioRel_max_over_AllComp), 1);
% 
%             % Evaluate the Posterior for each initial conditions
%             theta_Trials = theta_store_Mean_List_cell_k';
%             [Fields_Trials, GradxFields_Trials, GradyFields_Trials] = convert_theta_U1K0_to_Fields(theta_Trials);  % Should be an inline function
%             % Note: each row corresponds to a sample initial condition
%             
%             % Only use the first level
%             Nlevel = size(TrajJumps, 2);
%             TrajJumps_cellk = TrajJumps(cell_k, :);
%             logPost_cellk_Trials = zeros(size(theta_Trials, 1), 1);
%             
%             for nTrial_ind = 1:size(theta_Trials, 1)
%                 theta_cellk = theta_Trials(nTrial_ind, :);
%                 Fields_cellk = Fields_Trials(nTrial_ind, :);
%                 GradxFields_cellk = GradxFields_Trials(nTrial_ind, :);
%                 GradyFields_cellk = GradyFields_Trials(nTrial_ind, :);
%                 
%                 logLik_cellk = 0;
%                 for level = 1:Nlevel
%                     TrajJumps_cellk_level = TrajJumps_cellk(level);
%                     
%                     logLik_cellk = logLik_cellk + compute_logLik_jump_cell_linearSDE_vectorised(TrajJumps_cellk_level, Fields_cellk, GradxFields_cellk, GradyFields_cellk);
%                 end
%                 
%                 logPrior_cellk = compute_logPrior_cell_KPolar(Prior_Param, theta_cellk);
% 
%                 logPost_cellk_Trials(nTrial_ind) = logPrior_cellk + logLik_cellk;
%             end
%             
%             nTrial_tmp = find(logPost_cellk_Trials == max(logPost_cellk_Trials), 1);
%             
%             
%             if (nTrial_maxAccRat ~= nTrial_tmp)
%                   fprintf('----------------------------------------- \n')
% %                   theta_store_Mean_all = theta_store_Mean_List_cell_k'
% % %                 theta_SSD_Final_all = theta_SSD_List_cell_k'
% % %                 theta_MAP_tmp = theta_MAP(cell_k, :)
% %                 logPost_cellk_Trials_all = logPost_cellk_Trials'
%                 fprintf('Picked nTrial %d instead of %d in cell %d. \n', nTrial_tmp, nTrial_maxAccRat, cell_k)
%             else
%                 fprintf('----------------------------------------- \n')
%                 fprintf('----------------------------------------- \n')
% %                 theta_store_Mean_all = theta_store_Mean_List_cell_k'
% %                 % theta_SSD_Final_all = theta_SSD_List_cell_k'
% %                 theta_MAP_tmp = theta_MAP(cell_k, :)
% %                 logPost_cellk_Trials_all = logPost_cellk_Trials'
%                 fprintf('No change of picking nTrial %d and %d in cell %d. \n', nTrial_tmp, nTrial_maxAccRat, cell_k)
%             end
%             
%             % Set the theta_SSD for sigma1 and sigma2 to be equal
%             theta_SSD_Picked_cell_k =  theta_SSD_List_cell_k(:, nTrial_tmp);
%             theta_SSD_Picked_cell_k(7:8) = mean(theta_SSD_Picked_cell_k(7:8));
%             
%             theta_Picked(cell_k, :) = theta_store_Mean_List_cell_k(:, nTrial_tmp);
%             AccRatio_Pick(cell_k, :) = AccRatio_List_cell_k(:, nTrial_tmp);
%             theta_SSD_Picked(cell_k, :) = theta_SSD_Picked_cell_k;
%             nTrial_Picked(cell_k,:) = nTrial_tmp;
%     end
%     
% %   %%
% %   % START: Post-Process: Spatial smoothening with median filter
%     % 1a. For theta_init
% %     theta_Picked_ij = convert_cellk_to_cellij( theta_Picked, RectMesh_Param );
% % 
% %     theta_Picked_ij_FILTERED = median_filter_Dataijk(theta_Picked_ij);
% 
%     theta_Picked_FINAL = theta_Picked;
%     
% %     theta_Picked_k_FILTERED = convert_cellij_to_cellk(theta_Picked_ij_FILTERED, RectMesh_Param);
% %     theta_Picked_FINAL(:,1:2) = theta_Picked_k_FILTERED(:,1:2);
% 
%     % 2.  For theta_SSD
%     theta_SSD_Picked_ij = convert_cellk_to_cellij( theta_SSD_Picked, RectMesh_Param );
%     theta_SSD_Picked_ij_FILTERED = median_filter_Dataijk(theta_SSD_Picked_ij);
% 
%     theta_SSD_Picked_FINAL = convert_cellij_to_cellk(theta_SSD_Picked_ij_FILTERED, RectMesh_Param);
%     
% % %     % Set up theta_Init and MHRW_Param for actual inference
% % %     grid_vis = setup_grid_vis(Mesh, RectMesh_Param);
% % % 
% % %     figure(1)
% % %     ax_sp = plot_dataij_k(grid_vis, theta_Picked_ij);
% % %     figure(2)
% % %     ax_sp = plot_dataij_k(grid_vis, theta_Picked_ij_FINAL);
% % %     figure(11)
% % %     ax_sp = plot_dataij_k(grid_vis, log10(theta_SSD_Picked_ij));
% %     figure(12)
% %     ax_sp = plot_dataij_k(grid_vis, (theta_SSD_Picked_ij_Final));
% %     
% %     theta_Picked = convert_cellij_to_cellk( theta_Picked_ij_FINAL, RectMesh_Param );
% %     theta_SSD_Picked = convert_cellij_to_cellk( theta_Picked_ij_FINAL, RectMesh_Param );
% %     
% %     theta_Picked_FINAL_vis = convert_cellk_to_cellij(theta_Picked_FINAL, RectMesh_Param);
% %     ax_sp = plot_dataij_k(grid_vis, (theta_Picked_FINAL_vis));