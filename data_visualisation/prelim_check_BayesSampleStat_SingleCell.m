%% A quick plot visualisation of one cell inference 
close all
for expt = 1:2;

%theta_store_zip = theta_Stat_List{expt}.theta_store(:,:,[1,2:9:end]);
Ntheta_store = size(theta_Stat_List{expt}.theta_store, 3);
theta_store = theta_Stat_List{expt}.theta_store(:,:,round(0.5*Ntheta_store):end);
%theta_store = theta_store_zip;
[Fields, GradxFields, GradyFields] = convert_theta_to_Fields(theta_store);
trA = 0.5*(GradxFields(:,1,:) + GradyFields(:,2,:));

% MAP
theta_MAP = theta_Stat_List{expt}.theta_MAP;
[Fields_MAP, GradxFields_MAP, GradyFields_MAP] = convert_theta_to_Fields(theta_MAP);


%close all
Fields_mean = mean(Fields, 3);
GradxFields_mean = mean(GradxFields, 3);
GradyFields_mean = mean(GradyFields, 3);

A_mean = [[GradxFields_mean(1), GradyFields_mean(1)]; [GradxFields_mean(2), GradyFields_mean(2)]];
b_mean = [Fields_mean(1); Fields_mean(2)];

kappa_mean = [Fields_mean(3), Fields_mean(5); Fields_mean(5), Fields_mean(4)];


A_MAP = [[GradxFields_MAP(1), GradyFields_MAP(1)]; [GradxFields_MAP(2), GradyFields_MAP(2)]]
b_MAP = [Fields_MAP(1); Fields_MAP(2)]

kappa_MAP = [Fields_MAP(3), Fields_MAP(5); Fields_MAP(5), Fields_MAP(4)]
 
figure(1)
histogram(Fields(1,1,:))
hold on
histogram(Fields(1,2,:))

figure(2)
histogram(Fields(1,3,:))
hold on
histogram(Fields(1,4,:))
histogram(Fields(1,5,:))

figure(3)
histogram(GradxFields(1,1,:))
hold on
histogram(GradxFields(1,2,:))

figure(4)
histogram(GradyFields(1,1,:))
hold on
histogram(GradyFields(1,2,:))

figure(5)
histogram(trA(1,1,:))
hold on

figure(11)
histogram(theta_store(1,1,:))
hold on

figure(12)
histogram(theta_store(1,7,:))
hold on
histogram(theta_store(1,8,:))


figure(13)
histogram(theta_store(1,2,:))
hold on
histogram(theta_store(1,5,:))
histogram(theta_store(1,9,:))

figure(14)
histogram(theta_store(1,3,:))
hold on
histogram(theta_store(1,4,:))

figure(15)
histogram(theta_store(1,6,:))
hold on
end
