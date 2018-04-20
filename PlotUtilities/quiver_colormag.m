function qv_str = quiver_colormag(ax_cur, X, Y, U_vis, V_vis)
qv_str= quiver(ax_cur, X, Y, U_vis, V_vis);
jet_cmap = colormap('jet');

UV_mag = sqrt(U_vis.^2 + V_vis.^2);
%set(qv_str,'AutoScale','on', 'AutoScaleFactor',0.5)

%// Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(UV_mag, size(jet_cmap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), jet_cmap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(qv_str.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(qv_str.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
end