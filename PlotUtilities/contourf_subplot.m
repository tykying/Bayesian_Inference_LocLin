function ax = contourf_subplot(grid_vis, contourf_vis_data, subplot_param)
assert(all(size(contourf_vis_data) == size(subplot_param)))
assert(length(subplot_param(1).Position)==4);
% Default: [0 1 0 1]
% ax_han = axes('Units','Normalized','Position', subplot_param(1).Position);

fig_xI = subplot_param(1,1).Position(1);
fig_yI = subplot_param(1,1).Position(2);
fig_dx = subplot_param(1,1).Position(3);
fig_dy = subplot_param(1,1).Position(4);

% Mind the dimension: matrix size dimension is different from visualisation
nPlot_x = size(contourf_vis_data, 2);  % Not the first dimension!
nPlot_y = size(contourf_vis_data, 1);

% Set parameters
subplot_margin =  0.00625;

subplots_xI = fig_xI + subplot_margin;
subplots_xF = fig_xI + fig_dx - subplot_margin;
subplots_dx = (subplots_xF - subplots_xI - nPlot_x*subplot_margin*2)/nPlot_x;

subplots_yI = fig_yI + subplot_margin;
subplots_yF = fig_yI + fig_dy - subplot_margin;
subplots_dy = (subplots_yF - subplots_yI - nPlot_y*subplot_margin*2)/nPlot_y;
            
for j_x = 1:nPlot_x  % j <-> number of plots along the horizontal direction
    for i_y = 1:nPlot_y
        subplot_k_xI = fig_xI+(j_x-1)*(subplots_dx+2*subplot_margin) + subplot_margin;
        subplot_k_yI = fig_yI+((nPlot_y-i_y+1)-1)*(subplots_dy+2*subplot_margin) + subplot_margin;
        
        ax{i_y,j_x} = subplot('Position', [subplot_k_xI subplot_k_yI subplots_dx subplots_dy]);
        
        ContourLineValue_data = subplot_param(i_y, j_x).ContourLineValue_data;
        caxis_data = subplot_param(i_y, j_x).caxis_data;
        title_data = subplot_param(i_y, j_x).title_data;
        colormap_data = subplot_param(i_y, j_x).colormap_data;
        
        contourf(ax{i_y, j_x}, grid_vis.X, grid_vis.Y, contourf_vis_data{i_y, j_x}, ContourLineValue_data);
        caxis(ax{i_y, j_x}, caxis_data);
        title(ax{i_y, j_x}, title_data)

        pbaspect(ax{i_y,j_x}, [1 1 1])
        set(ax{i_y,j_x}, 'xticklabel', []);
        set(ax{i_y,j_x}, 'yticklabel', []);
        colormap(ax{i_y,j_x}, colormap_data);
        colorbar(ax{i_y,j_x}, 'Location', 'southoutside');        
    end
end

% Older version; contains more special treatments
%assert(all(size(contourf_vis_data) == size(title_vis_data)))
% NC_Subplot = numel(contourf_vis_data);
% subplot_iterator = 1;
% for row_ind = 1:size(contourf_vis_data, 1)
%     for col_ind = 1:size(contourf_vis_data, 2)
%         %subplot_string = ['Sampling Interval = ', num2str(SamplingInterval_List(iterator)/time_scale), '; '];
%         
%         X_vis = grid_vis.X;
%         Y_vis = grid_vis.Y;
%         
%         subplot(size(contourf_vis_data, 1), size(contourf_vis_data, 2), subplot_iterator)
%         
%         colormap(jet);
%         
%         % Division by Sampling Interval
%         Field_vis = contourf_vis_data{row_ind, col_ind};
%         
%         % Set ContourLineValue
%         ContourLineValue = linspace(min(Field_vis(:)), max(Field_vis(:)), 15);
% 
%         if strcmp(subplot_param.caxis_data{row_ind, col_ind}, 'sym')
%             Field_vis_AbsMax = max(abs(Field_vis(:)));
%             ContourLineValue = linspace(-Field_vis_AbsMax, Field_vis_AbsMax, 10);
%         end
%         if ~isempty(subplot_param.ContourLineValue_data{row_ind, col_ind})
%             ContourLineValue = subplot_param.ContourLineValue_data{row_ind, col_ind};
%         end
%         
%         contourf(X_vis, Y_vis, Field_vis, ContourLineValue);
%         hold on
% 
%         % Set caxis        
%         if ~isempty(subplot_param.caxis_data{row_ind, col_ind})
%             if strcmp(subplot_param.caxis_data{row_ind, col_ind}, 'sym')
%                 caxis([-Field_vis_AbsMax, Field_vis_AbsMax])
%             elseif size(subplot_param.caxis_data{row_ind, col_ind} == 2)
%                 caxis(subplot_param.caxis_data{row_ind, col_ind})
%             end
%             
%         end
%        
%         
%         % Set Aspect Ratio
%         if isempty(subplot_param.set_pbaspect{row_ind, col_ind})
%             pbaspect([1 1 1])
%             set(gca, 'xticklabel', []);
%             set(gca, 'yticklabel', []);
%         end
%                 
%         title(subplot_param.title_vis_data{row_ind, col_ind});
%         subplot_iterator = subplot_iterator + 1;
%                 
%         if (isempty(subplot_param.display_globalCbar{1}))
%             c = colorbar;
%             c.Location = 'southoutside';
%         end
%     end
% end
% 
% if (subplot_param.display_globalCbar{1} == 1)
%     c = colorbar;
%     c.Location = 'southoutside';
%     c.Position = [0.25 0.05 0.5 0.02];
% end

end