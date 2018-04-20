function [plot_obj, text_obj] = plot_MarkCellCentre(grid_vis, cell_ij_list)
    % Mark the points
    DX = grid_vis.X(2,1) - grid_vis.X(1,1);
    DY = grid_vis.Y(1,2) - grid_vis.Y(1,1);
    
    for list_iterator = 1:length(cell_ij_list(:))
        cell_i = cell_ij_list{list_iterator}(1);
        cell_j = cell_ij_list{list_iterator}(2);
        
        x_coord = grid_vis.cell_centre(cell_i, cell_j, 1);
        y_coord = grid_vis.cell_centre(cell_i, cell_j, 2);
        
        %plot(x,y,'d')
        plot_obj = plot(x_coord, y_coord, 'Color','k', 'Marker','square', 'MarkerSize', 16);
        %text_obj = text(x_coord + DX*(0.15), y_coord - DY*(0.25), sprintf('(%i, %i)', cell_i, cell_j), 'Color','k', 'FontSize', 18);
        text_obj = text(x_coord - 2*DX, y_coord + DY, sprintf('(%i, %i)', cell_i, cell_j), 'Color',0.2*[1,1,1], 'FontSize', 14);

    end
end