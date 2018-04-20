function qv_str = quiver_angle(ax_cur, X_vis, Y_vis, U_vis, V_vis)
UV_mag = sqrt(U_vis.^2 + V_vis.^2);

cos_Unorm = U_vis./UV_mag;
sin_Vnorm = V_vis./UV_mag;
white_vector = [1,1,1];

Dx = unique(diff(X_vis, 1, 1));
Dy = unique(diff(Y_vis, 1, 2));

% Ensure square grid
assert(length(Dx) == length(Dy));
assert(length(Dx) == 1);

% Shift the starting point of the arrows
X_vis_lv = X_vis - 0.5*cos_Unorm.*Dx;
Y_vis_lv = Y_vis - 0.5*sin_Vnorm.*Dy;

qv_str= quiver(ax_cur, X_vis_lv, Y_vis_lv, U_vis, V_vis, 'color',0.4*white_vector);
set(qv_str,'AutoScale','on', 'AutoScaleFactor',0.75); %, 'AlignVertexCenters', 'on')
end