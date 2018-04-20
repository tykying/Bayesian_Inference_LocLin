function [ax_quiver] = visualise_LinearVelocityFld(theta_cell_k, bin_sizes, SamplingInterval_vis)
%%
A = [theta_cell_k(3), theta_cell_k(4); theta_cell_k(5), theta_cell_k(6)];
b = [theta_cell_k(1); theta_cell_k(2)];

x_GP = linspace(-0.5*bin_sizes(1), 0.5*bin_sizes(1), 8+1);
y_GP = linspace(-0.5*bin_sizes(2), 0.5*bin_sizes(2), 8+1);

[x,y] = meshgrid(x_GP,y_GP);

u = zeros(size(x));
v = zeros(size(x));
for k = 1:length(x(:))
    X_k = [x(k); y(k)];
    
    U_k = A*X_k + b;
    
    u(k) = U_k(1);
    v(k) = U_k(2);
end

ax_quiver = quiver(x,y,u,v);

% Release a particle at cell centre
h = SamplingInterval_vis*3600*24;
%h = 128*3600*24;
tspan = [0 h];

%Y0 = [];
SeedPointInd = linspace(-1, 1, 64+1);
for i = SeedPointInd
    for j = SeedPointInd
        %Y0 = [Y0 [i; j] .* 0.5.* bin_sizes'];
        Y0 = [i; j] .* 0.5.* bin_sizes';
        [t,y] = ode45(@(t,x) A*x + b, tspan, Y0);
        
        hold on
        %scatter(y([1,end], 1), y([1,end], 2))
        plot(y(:, 1), y(:, 2), 'g');
    end
end