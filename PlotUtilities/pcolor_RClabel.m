function ax_sp = pcolor_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List)
ax_sp = cell(size(Data_ij_vis_List));
Nrow = size(Data_ij_vis_List, 1);
Ncol = size(Data_ij_vis_List, 2);

assert(Nrow == length(YLabelString_List));
assert(Ncol == length(TitleString_List));

% Initialise ax_sp
ind = 1;
for r_ind = 1: Nrow
    for c_ind = 1:Ncol
        ax_sp{r_ind, c_ind} = subplot(size(ax_sp, 1), size(ax_sp, 2), ind);
        ind = ind + 1;
    end
end

for c_ind = 1:Ncol
    for r_ind = 1:Nrow
        ax_spc = ax_sp{r_ind, c_ind};
        surface_sp = pcolor(ax_spc, grid_vis.X, grid_vis.Y , Data_ij_vis_List{r_ind, c_ind} );  
        %surface_sp = contourf(ax_spc, grid_vis.X, grid_vis.Y , Data_ij_vis_List{r_ind, c_ind} );
        %set(ax_spc,'xticklabel',{[]}); set(ax_spc,'yticklabel',{[]});
        if isempty(row_CLim{r_ind}) == 0
            ax_spc.CLim = row_CLim{r_ind};
        end
        %colorbar(ax_spc);
        %pbaspect(ax_spc, [1 1 1]);
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

for ind = 1:length(ax_sp(:))
    subplotax = subplot(Nrow, Ncol, ind);
    run('Script_AxesConfig.m')
end

end