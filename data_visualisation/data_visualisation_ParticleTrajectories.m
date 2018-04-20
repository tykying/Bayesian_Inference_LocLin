%% Plot Particle Trajectories

space_scale = 1000;

% Given x, y, ts_list in the memory
veloc_testcase = 'TaylorGreen';
TrajRawDataPath = '/home/s1046972/opt/qgm2_particle/PART_TRAJ/traj_Npart2000_tIntv9yr_SampletIntv24hr.mat'
load(TrajRawDataPath)

%%% Convert x, y from centimeters to meters
x = x/(100*space_scale); y = y/(100*space_scale);  L = L/(100*space_scale);
h = (ts_list(2)-ts_list(1));
    
space_scale = 1000;
time_scale = 30*24*3600;

veloc_testcase = 'TaylorGreen';
TrajRawDataPath = '/home/s1046972/Desktop/Two_Dimensional/Trajectories_Generated/taylor_green_trajectories_saved.mat'
load(TrajRawDataPath)

SamplingInterval = 0.1;
SamplingParticle = 40;

SamplingParticle = 10;
time_scale = 1;


[Jumps, TrajData_SubSampled] = convert_traj_to_JumpsStruct(x, y, ts_list, SamplingInterval, SamplingParticle);

video = VideoWriter(['./Video_Generated/', veloc_testcase, '_traj.avi'],'Uncompressed AVI');
video.FrameRate = 10;
open(video)

f11 = figure(11);
set(f11, 'Position', [100, 100, 800, 800]);
for i = 1:length(TrajData_SubSampled.ts_list_SubSampled(1:1000))
tn = TrajData_SubSampled.ts_list_SubSampled(i) - TrajData_SubSampled.ts_list_SubSampled(1);
c = linspace(1,25,length(TrajData_SubSampled.x_SubSampled(i,:)));
%plot_MeanFlowStreamline(space_scale);

scatter(TrajData_SubSampled.x_SubSampled(i,:), TrajData_SubSampled.y_SubSampled(i,:), [], c, 'filled');

tn_string = sprintf('%4.2f', tn/time_scale);

xlim([-30, 30])
ylim([-30, 30])
%title(['time = ', tn_string, ' months'])
%xlabel('x (in km)');
%xlabel('y (in km)');
title(['t = ', tn_string])
xlabel('x');
xlabel('y');

drawnow();
frame = getframe(f11);
writeVideo(video,frame);

clf
end
close(video)


f12 = figure(12);
set(f12, 'Position', [100, 100, 800, 800]);
plot(TrajData_SubSampled.x_SubSampled(1:1000, :), TrajData_SubSampled.y_SubSampled(1:1000, :))


xlim([0, 4000])
ylim([0, 4000])
xlabel('x (in km)');
ylabel('y (in km)');

drawnow();

clf
