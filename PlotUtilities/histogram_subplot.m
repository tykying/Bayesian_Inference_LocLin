function histogram_subplot(histogram_vis_data, title_vis_data)

assert(all(size(histogram_vis_data) == size(title_vis_data)))
NC_Subplot = numel(histogram_vis_data);

subplot_iterator = 1;
for row_ind = 1:size(histogram_vis_data, 1)
    for col_ind = 1:size(histogram_vis_data, 2)
        %subplot_string = ['Sampling Interval = ', num2str(SamplingInterval_List(iterator)/time_scale), '; '];
        
        subplot(size(histogram_vis_data, 1), size(histogram_vis_data, 2), subplot_iterator)
                
        % Division by Sampling Interval
        Data_vis = histogram_vis_data{row_ind, col_ind};
        
        Hist_vis = histogram(Data_vis)
        Hist_vis.Normalization = 'pdf';    
        
        title(title_vis_data{row_ind, col_ind});
        subplot_iterator = subplot_iterator + 1;
        
        hold on
    end
end

end