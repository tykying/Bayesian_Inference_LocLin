function plotxy_subplot(plotx_vis_data, ploty_vis_data, title_vis_data)

assert(all(size(plotx_vis_data) == size(title_vis_data)))
assert(all(size(ploty_vis_data) == size(title_vis_data)))


subplot_iterator = 1;
for row_ind = 1:size(title_vis_data, 1)
    for col_ind = 1:size(title_vis_data, 2)
        %subplot_string = ['Sampling Interval = ', num2str(SamplingInterval_List(iterator)/time_scale), '; '];
        
        subplot(size(title_vis_data, 1), size(title_vis_data, 2), subplot_iterator)
                
        % Division by Sampling Interval
        plotx_vis = plotx_vis_data{row_ind, col_ind};
        ploty_vis = ploty_vis_data{row_ind, col_ind};

        plot(plotx_vis, ploty_vis)
        
        title(title_vis_data{row_ind, col_ind});
        subplot_iterator = subplot_iterator + 1;
        
        % ad-hoc: should read the time scale
        %xlabel('Days')
        %ylabel('Autocorrelation')
        hold on
    end
end

end