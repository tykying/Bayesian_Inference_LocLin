function [imagesc_sp, ax_sp] = plot_dataij_k(grid_vis, dataij_k)
    NPlot = size(dataij_k, 3);
    NRow = fix(NPlot/3-10E-9)+1;
    
    Title_List = {'u', 'v', 'd_x u', 'd_y u', 'd_x v', 'd_y v', '\sigma_1', '\sigma_2', '\phi'};
    
    for k = 1:NPlot
        ax_sp(k) = subplot(NRow,3,k);
        dataij = dataij_k(:, :, k);
        %pcolor(grid_vis.X, grid_vis.Y, dataij);
        
        data05 = prctile(dataij(:),5);
        data95 = prctile(dataij(:),95);
        
        
        %imagesc_sp = imagesc(ax_sp(k), grid_vis.X(:,1), grid_vis.Y(1,:), dataij');
        X_Ind = 1:length(grid_vis.X(:,1));
        Y_Ind = 1:length(grid_vis.Y(1,:));
        imagesc_sp(k) = imagesc(ax_sp(k), X_Ind, Y_Ind, dataij','AlphaData',~isnan(dataij'));
        
        
        set(ax_sp(k),'xticklabel',{[]});
        set(ax_sp(k),'yticklabel',{[]});
        axis(ax_sp(k), 'xy')
        if NPlot == 9
            title(Title_List{k})
        end
            
        %data_ub = 4.0*median(abs(dataij(:)));
        caxis([data05, data95]);
        %title(sprintf('Cmp %d', k));
        colormap('jet')
        colorbar
        axis square
    end
end