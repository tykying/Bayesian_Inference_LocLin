function ax_sp = plot_RClabel(x, Data_ij_vis_List, YLabelString_List, TitleString_List)
ax_sp = cell(size(Data_ij_vis_List));
Nrow = size(Data_ij_vis_List, 1);
Ncol = size(Data_ij_vis_List, 2);

assert(Nrow == length(YLabelString_List));
assert(Ncol == length(TitleString_List));

MaxData = zeros(size(Data_ij_vis_List));
MinData = zeros(size(Data_ij_vis_List));

% Initialise ax_sp
ind = 1;
for r_ind = 1: Nrow
    for c_ind = 1:Ncol
        ax_sp{r_ind, c_ind} = subplot(size(ax_sp, 1), size(ax_sp, 2), ind);
        ind = ind + 1;
    end
end

ylim_max = zeros(1,2);
for c_ind = 1:Ncol
    for r_ind = 1:Nrow
        ax_spc = ax_sp{r_ind, c_ind};
        plot_sp = plot( ax_spc, x, Data_ij_vis_List{r_ind, c_ind} );
        
        MaxData(r_ind, c_ind) = max(Data_ij_vis_List{r_ind, c_ind}(:));
        MinData(r_ind, c_ind) = min(Data_ij_vis_List{r_ind, c_ind}(:));
        
        %set(ax_spc,'xticklabel',{[]}); set(ax_spc,'yticklabel',{[]});
        %         if isempty(row_CLim{r_ind}) == 0
        %             ax_spc.CLim = row_CLim{r_ind};
        %         end
        %colorbar(ax_spc);
        %pbaspect(ax_spc, [1 1 1]);
        set(ax_spc, 'XLim', x([1,end]));
        
        % Set identical ylim for all figures
        % Obtain ylim
        ylim_max(1) = min(ax_spc.YLim(1), ylim_max(1));
        ylim_max(2) = max(ax_spc.YLim(2), ylim_max(2));
    end
end


% Labelling: Along each column
for c_ind = 1:Ncol
    ax_spc = ax_sp{1, c_ind};
    ax_spc.Title.String = TitleString_List{c_ind};
end

% Labelling: Along each row
for r_ind = 1:Nrow
    ax_sp{r_ind, 1}.YLabel.String = YLabelString_List{r_ind};
    
    subplotax = subplot(Nrow, Ncol, Ncol*r_ind);
    ax_sp_pos = get(subplotax,'Position');
    %colorbar('Position', [ax_sp_pos(1)+ax_sp_pos(3)+0.01  ax_sp_pos(2)  0.01  ax_sp_pos(4)]);
end

% Amend axis
globalylim = 0;  % Use the same ylim for all subplots
for c_ind = 1:Ncol
    for r_ind = 1:Nrow
        ax_spc = ax_sp{r_ind, c_ind};       
        
        if globalylim == 1
            ylim_vis = ylim_max;
        else
            MaxData_ind = MaxData(r_ind, c_ind);
            MinData_ind = MinData(r_ind, c_ind);
            ylim_vis = obtain_ylim_OrdMag(MaxData_ind, MinData_ind);
        end
        
        set(ax_spc, 'YLim', 1.05*ylim_vis);
    end
end

for ind = 1:length(ax_sp(:))
    subplotax = subplot(Nrow, Ncol, ind);
    run('Script_AxesConfig.m');
end

end