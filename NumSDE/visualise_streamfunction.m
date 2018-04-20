%% Below is some ad-hoc code to visualise streamfunction
% k=4;
% npart = 32;
% pid = linspace(1, 1024, npart)
% plot(trajdata_Pe(k).x(:, pid), trajdata_Pe(k).y(:, pid));
% title(num2str(Pe_List(k)));
% ylim([-10, 10])
% xlim([-500, 500])

L_D = 8;

xgv = linspace(-L_D, L_D, 256);
ygv = linspace(-L_D, L_D, 256);

[X,Y] = meshgrid(xgv,ygv);
grid_vis = struct('X', X, 'Y', Y);

psi_scale = 0.5;
Lx = 1;
Ly = 1;
Ls = 8/(2*pi);
%Ls = 4;

phi = pi/4*1;
scale = 0.4;
b1 = scale*cos(phi);
b2 = scale*sin(phi);

psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y) + 0.5*(-Au1(2)*y.^2 + Au2(1)*x.^2) -2*Au1(1)*x.*y;
psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y) - (4*Ls)*sin(y/Ls);
psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y) - b1*y + b2*x;


% Disable for video production
%PSI = psi(X, Y);
%contourf(X, Y, PSI)

% hold on
% psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y);
% 
% PSI = psi(X, Y);
% contour(X, Y, PSI)


%% 

% Make video of the P.D.F. of Jumps at different sampling interval
TrajType_Save = 'Streamfunction';
veloc_testcase = 'tg_w_mean';

u_scale_List = linspace(0, 1.28, 128+1);  % Ad-hoc
u_angle_List = linspace(0, pi, 128+1);

%u_scale_List = 0.64;  % Ad-hoc
u_angle_List = pi/4;


vsU = 1;
if vsU == 1
    u_angle_List = [pi/4*0.5];  % Angle of incident mean flow
    veloc_testcase = [veloc_testcase, '_phi', num2str(u_angle_List)];
else
    u_scale_List = [0.16];
    veloc_testcase = [veloc_testcase, '_U', num2str(u_scale_List)];
end


video = VideoWriter(['/data/tying/Video_Generated/', TrajType_Save, '_', veloc_testcase, '_vsU.avi'],'Uncompressed AVI');
video.FrameRate = 10
open(video)
vframe = figure(22)
colormap(jet(256))
set(vframe, 'Position', [100, 100, 800, 800]);
for phi = u_angle_List
for u_scale = u_scale_List
    clf;

    b1 = u_scale*cos(phi);
    b2 = u_scale*sin(phi);
    
    psi = @(x, y) psi_scale*sin(Lx*x).*sin(Ly*y) - b1*y + b2*x;
    
    PSI = psi(X, Y);
    PSI = PSI/max(abs(PSI(:)));
    
    % Streamfunctions
    ContourLineValue = linspace(-1, 1, 16);
    contourf(X, Y, PSI, ContourLineValue)
    caxis([ContourLineValue(1), ContourLineValue(end)]);
    xlim([min(grid_vis.X(:)), max(grid_vis.X(:))]);
    ylim([min(grid_vis.Y(:)), max(grid_vis.Y(:))]);
    
    title(['U= ', num2str(u_scale), '; Angle=', num2str(phi)]);
    %colorbar;
    pbaspect([1 1 1])
    
    drawnow();
    frame = getframe(vframe);
    writeVideo(video,frame);
    
end
end
close(video)