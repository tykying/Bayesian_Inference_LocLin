FigObj = figure('Name','Velocities','NumberTitle','off', 'Units', 'normalized', 'pos', [0.2 0.2 0.6 0.6]);

contourf_vis_data = {u_vis, dudx_vis, dudy_vis; ...
    v_vis, dvdx_vis, dvdy_vis; };

if ~exist('U_CLV')
    U_CLV = []; dUdx_CLV = []; dUdy_CLV = [];
    V_CLV = []; dVdx_CLV = []; dVdy_CLV = [];
end
if ~exist('U_range')
    U_range = [-1,1]*max(abs(u_vis(:)));
    dUdx_range = [-1,1]*max(abs(dudx_vis(:)));
    dUdy_range = [-1,1]*max(abs(dudy_vis(:)));
    V_range = [-1,1]*max(abs(v_vis(:)));
    dVdx_range = [-1,1]*max(abs(dvdx_vis(:)));
    % dVdy_range = [-1,1]*max(abs(dvdy_vis(:)));
    dVdy_range = dUdx_range;    % Check Incompressiblity
end

ContourLineValue_data = {U_CLV, dUdx_CLV, dUdy_CLV; V_CLV, dVdx_CLV, dVdy_CLV};
caxis_data = {U_range, dUdx_range, dUdy_range; V_range, dVdx_range, dVdy_range};
colormap_data = {'jet', 'jet', 'jet'; 'jet', 'jet', 'jet'};

title_data = {'u', 'dudx', 'dudy'; 'v', 'dvdx', 'dvdy'};
subplot_param = struct('Position', figPos, ...
    'ContourLineValue_data', ContourLineValue_data, ...
    'caxis_data', caxis_data, 'title_data', title_data, ...
    'colormap_data', colormap_data);

pcolor_subplot(grid_vis, contourf_vis_data, subplot_param);